"""Single-channel natural-weighted test images per field per VLA configuration.

For every field and every VLA configuration present in the MS, run tclean
(niter=0) on the middle channel of the first SPW, selecting only that
configuration's observation IDs. Only the PSF and residual are kept.

Configurations are identified per OBSERVATION_ID from the 95th-percentile uv
distance compared with absolute nominal A/B/C/D baselines
(config_from_p95), not by rank.

Run with:
    casa --nologger --nogui --log2term -c run_test_imaging.py <ms> <outdir>

Writes <outdir>/test_<field>_<config>.{psf,residual} and
<outdir>/test_imaging_summary.txt.
"""

import os
import shutil
import sys

import numpy as np

from casatasks import tclean
from casatools import image as iatool
from casatools import msmetadata, table


def script_args(script_name):
    """Arguments following `script_name` on a `casa -c script args` command line."""
    idx = max(i for i, a in enumerate(sys.argv) if a.endswith(script_name))
    return sys.argv[idx + 1:]


# ---------------------------------------------------------------------------
# VLA configuration identification (absolute, not rank-based)
# ---------------------------------------------------------------------------

# Nominal maximum baselines (m). The 95th-percentile uv distance of a track is
# ~0.7 x the maximum baseline, so compare p95 against 0.7 x nominal.
NOMINAL_MAX_BASELINE = {"A": 36.4e3, "B": 11.1e3, "C": 3.4e3, "D": 1.03e3}
P95_FRACTION = 0.7


def config_from_p95(p95):
    """Map a p95 uv distance (m) to A/B/C/D by the log-midpoint between configs."""
    order = ["A", "B", "C", "D"]
    nominal = {k: P95_FRACTION * NOMINAL_MAX_BASELINE[k] for k in order}
    for hi, lo in zip(order[:-1], order[1:]):
        if p95 >= np.sqrt(nominal[hi] * nominal[lo]):
            return hi
    return "D"


def obsid_p95_baselines(tb, vis, chunk_size=500000, uvmin=10.0, uvmax=1.0e5, nbins=400):
    """Per-OBSERVATION_ID p95 uv distance (m) from unflagged rows (FLAG_ROW only).

    `tb` is the CASA table tool. Uses a log-binned histogram so memory is bounded.
    """
    edges = np.logspace(np.log10(uvmin), np.log10(uvmax), nbins + 1)
    hist = {}
    tb.open(vis)
    nrows = tb.nrows()
    for start in range(0, nrows, chunk_size):
        n = min(chunk_size, nrows - start)
        obsid = tb.getcol("OBSERVATION_ID", startrow=start, nrow=n)
        uvw = tb.getcol("UVW", startrow=start, nrow=n)
        frow = tb.getcol("FLAG_ROW", startrow=start, nrow=n)
        uvd = np.hypot(uvw[0], uvw[1])
        good = (~frow) & (uvd > 0)
        idx = np.clip(np.digitize(uvd, edges) - 1, 0, nbins - 1)
        for oid in np.unique(obsid[good]):
            sel = good & (obsid == oid)
            h = hist.setdefault(int(oid), np.zeros(nbins, dtype=np.int64))
            h += np.bincount(idx[sel], minlength=nbins)
    tb.close()

    out = {}
    for oid, h in hist.items():
        cum = np.cumsum(h) / h.sum()
        i = min(int(np.searchsorted(cum, 0.95)), nbins - 1)
        lo = cum[i - 1] if i > 0 else 0.0
        frac = 0.0 if cum[i] == lo else (0.95 - lo) / (cum[i] - lo)
        out[oid] = edges[i] + frac * (edges[i + 1] - edges[i])
    return out


def identify_configs(tb, vis):
    """Return ({config: [obsids]}, {config: p95 baseline (m)}, {obsid: p95})."""
    p95 = obsid_p95_baselines(tb, vis)
    groups = {}
    for oid, v in sorted(p95.items()):
        groups.setdefault(config_from_p95(v), []).append(oid)
    group_p95 = {c: float(np.max([p95[o] for o in oids])) for c, oids in groups.items()}
    return groups, group_p95, p95


# ---------------------------------------------------------------------------
# Imaging-size helper
# ---------------------------------------------------------------------------

def good_imsize(n, max_size=4096):
    """Smallest even number >= n whose only prime factors are 2, 3, 5 (capped)."""
    n = int(min(max(n, 64), max_size))
    while n <= max_size:
        m = n
        for f in (2, 3, 5):
            while m % f == 0:
                m //= f
        if m == 1 and n % 2 == 0:
            return n
        n += 1
    return max_size


C_LIGHT = 299792458.0
primary_beam_arcmin = 30.0   # field of view to cover (HI at 1.4 GHz is ~30')
max_imsize = 4096
cell_oversample = 5.0        # pixels per (lambda / p95 baseline)

args = script_args("run_test_imaging.py")
vis = args[0].rstrip("/")
outdir = args[1] if len(args) > 1 else os.path.dirname(os.path.abspath(vis)) + "/test_images"
os.makedirs(outdir, exist_ok=True)

tbl = table()
groups, group_p95, p95 = identify_configs(tbl, vis)
print("Config groups (obsids):", groups)
print("p95 baselines per obsid (m):", {k: round(v) for k, v in p95.items()})

tbl.open(vis)
datacolumn = "corrected" if "CORRECTED_DATA" in tbl.colnames() else "data"
tbl.close()

tbl.open(vis + "/DATA_DESCRIPTION")
spw = int(sorted(set(tbl.getcol("SPECTRAL_WINDOW_ID")))[0])
tbl.close()

msmd = msmetadata()
msmd.open(vis)
fields = list(msmd.fieldnames())
nchan = int(msmd.nchan(spw))
# chan = nchan // 2
chan = 0
freq_hz = float(msmd.chanfreqs(spw)[chan])
# Only image fields that have scans (skip unobserved/empty names)
fields_with_data = [f for f in dict.fromkeys(fields) if len(msmd.scansforfield(f)) > 0]
msmd.close()
wavelength = C_LIGHT / freq_hz
print(f"SPW {spw}: channel {chan} of {nchan}, {freq_hz / 1e6:.3f} MHz, datacolumn={datacolumn}")

rows = []
for config, obsids in sorted(groups.items()):
    baseline = group_p95[config]
    cell_arcsec = np.degrees(wavelength / baseline / cell_oversample) * 3600.0
    imsize = good_imsize(primary_beam_arcmin * 60.0 / cell_arcsec, max_imsize)
    obs_sel = ",".join(str(o) for o in obsids)
    for field in fields_with_data:
        prefix = os.path.join(outdir, f"test_{field}_{config}")
        row = dict(field=field, config=config, obsids=obs_sel, cell=cell_arcsec,
                   imsize=imsize, chan=chan, status="ok", psf_peak=np.nan, resid_rms=np.nan)
        try:
            for ext in (".psf", ".residual"):
                shutil.rmtree(prefix + ext, ignore_errors=True)
            tclean(vis=vis, imagename=prefix, field=field, observation=obs_sel,
                   spw=f"{spw}:{chan}", datacolumn=datacolumn, specmode="cube",
                   nchan=1, imsize=[imsize, imsize], cell=f"{cell_arcsec:.4f}arcsec",
                   gridder="standard", weighting="natural", niter=0,
                   deconvolver="hogbom", pblimit=0.7, restoration=False,
                   interactive=False, parallel=False)
            ia = iatool()
            ia.open(prefix + ".residual")
            row["resid_rms"] = float(ia.statistics()["rms"][0])
            ia.close()
            ia.open(prefix + ".psf")
            row["psf_peak"] = float(ia.statistics()["max"][0])
            ia.close()
        except Exception as exc:   # one bad field/config must not stop the rest
            row["status"] = f"FAILED: {exc}"
            print(f"tclean failed for {field} {config}: {exc}")
        finally:
            # Keep only the PSF and residual
            for fn in os.listdir(outdir):
                full = os.path.join(outdir, fn)
                if fn.startswith(os.path.basename(prefix) + ".") and \
                        not fn.endswith((".psf", ".residual")):
                    shutil.rmtree(full) if os.path.isdir(full) else os.remove(full)
        rows.append(row)

summary = os.path.join(outdir, "test_imaging_summary.txt")
with open(summary, "w") as fh:
    fh.write(f"Test imaging for {vis}\n")
    fh.write(f"SPW {spw} channel {chan}/{nchan} ({freq_hz / 1e6:.3f} MHz), natural weighting, "
             f"datacolumn={datacolumn}\n\n")
    fh.write("Config  obsids  p95_baseline_m\n")
    for c, o in sorted(groups.items()):
        fh.write(f"{c:>6}  {','.join(map(str, o)):>6}  {group_p95[c]:.0f}\n")
    fh.write("\n")
    fh.write("{:<16}{:>4}{:>12}{:>7}{:>10}{:>12}  status\n".format(
        "field", "cfg", "cell_arcsec", "imsize", "psf_peak", "resid_rms"))
    for r in rows:
        fh.write(f"{r['field']:<16}{r['config']:>4}{r['cell']:>12.3f}{r['imsize']:>7}"
                 f"{r['psf_peak']:>10.3f}{r['resid_rms']:>12.4g}  {r['status']}\n")
print(f"Wrote {summary}")
