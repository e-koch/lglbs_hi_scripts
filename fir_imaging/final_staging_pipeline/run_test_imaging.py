"""Single-channel natural-weighted test images per field per VLA configuration.

For every field and every VLA configuration present in the MS, run tclean
(niter=0) on channel 1 of ALL spectral windows in the MS (they are assumed to be the
same band on the same channel grid, e.g. the same line observed in different
correlator setups/epochs; checked below), combined into one output channel,
selecting only that configuration's observation IDs. Only the PSF and residual
are kept.

Configurations are identified per OBSERVATION_ID from the 95th-percentile uv
distance compared with absolute nominal A/B/C/D baselines
(config_from_p95), not by rank.

Run with:
    casa --nologger --nogui --log2term -c run_test_imaging.py <ms> <outdir> [mode [arg]]

Modes (the (field, config) jobs are independent, so they can be sharded over several
CASA processes; job1_prep_and_test.sh does this with TEST_SHARDS):
    (none)       plan + run every job serially + write the summary (original behaviour)
    plan         identify configs/fields once, write <outdir>/.test_imaging_plan.json
    shard I/N    run jobs k with k % N == I from the plan (I = 0..N-1), write
                 <outdir>/.test_imaging_shard_I.json. Needs only the plan, not the MS metadata.
    merge N      combine the N shard files into test_imaging_summary.txt (jobs with no
                 result, e.g. from a crashed worker, are listed as FAILED) and delete the
                 plan/shard files.

Writes <outdir>/test_<field>_<config>.{psf,residual} and
<outdir>/test_imaging_summary.txt (same format in every mode).
"""

import json
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

# The channel is chosen with tclean's start=<channel>, width=1 (no spw selector): an spw='0:1'
# selection could give blank PSFs (nothing gridded) for some MSs but not others.
# Channel 1, not 0: after combine_spws regridding, channel 0 can be (nearly) fully flagged
# for most fields (OH: ~98% of rows), which gave empty images for every mosaic pointing.
chan = 1

PLAN_NAME = ".test_imaging_plan.json"


def shard_name(i):
    return f".test_imaging_shard_{i}.json"


# ---------------------------------------------------------------------------
# Plan / run / summary
# ---------------------------------------------------------------------------

def build_plan(vis):
    """Identify configs, fields, SPWs and the per-(config, field) tclean settings."""
    tbl = table()
    groups, group_p95, p95 = identify_configs(tbl, vis)
    print("Config groups (obsids):", groups)
    print("p95 baselines per obsid (m):", {k: round(v) for k, v in p95.items()})

    tbl.open(vis)
    datacolumn = "corrected" if "CORRECTED_DATA" in tbl.colnames() else "data"
    tbl.close()

    tbl.open(vis + "/DATA_DESCRIPTION")
    spws = sorted(int(s) for s in set(tbl.getcol("SPECTRAL_WINDOW_ID")))
    tbl.close()

    msmd = msmetadata()
    msmd.open(vis)
    fields = list(msmd.fieldnames())
    nchans = {s: int(msmd.nchan(s)) for s in spws}
    nchan = min(nchans.values())
    chan_freqs = {s: float(msmd.chanfreqs(s)[chan]) for s in spws}
    chan_widths = {s: abs(float(msmd.chanwidths(s)[chan])) for s in spws}
    freq_hz = chan_freqs[spws[0]]
    # Only image fields that have scans (skip unobserved/empty names)
    fields_with_data = [f for f in dict.fromkeys(fields) if len(msmd.scansforfield(f)) > 0]
    msmd.close()
    wavelength = C_LIGHT / freq_hz
    spw_label = ",".join(str(s) for s in spws)
    print(f"SPW(s) {spw_label}: channel {chan} of {nchan}, {freq_hz / 1e6:.3f} MHz, datacolumn={datacolumn}")

    # Combining the SPWs into one channel assumes they share a channel grid. Warn loudly if not.
    spw_warning = None
    if len(spws) > 1:
        max_df = max(chan_freqs.values()) - min(chan_freqs.values())
        if len(set(nchans.values())) > 1 or max_df > min(chan_widths.values()) or \
                (max(chan_widths.values()) - min(chan_widths.values())) > 1e-3 * min(chan_widths.values()):
            spw_warning = (f"WARNING: SPWs {spw_label} do not share a channel grid "
                           f"(nchan={nchans}, ch{chan} spread {max_df:.0f} Hz, widths={chan_widths})")
            print(spw_warning)

    jobs = []
    for config, obsids in sorted(groups.items()):
        baseline = group_p95[config]
        cell_arcsec = float(np.degrees(wavelength / baseline / cell_oversample) * 3600.0)
        imsize = good_imsize(primary_beam_arcmin * 60.0 / cell_arcsec, max_imsize)
        obs_sel = ",".join(str(o) for o in obsids)
        for field in fields_with_data:
            jobs.append(dict(field=field, config=config, obsids=obs_sel, cell=cell_arcsec,
                             imsize=int(imsize), chan=chan))
    return dict(vis=vis, datacolumn=datacolumn, spw_label=spw_label, chan=chan, nchan=nchan,
                freq_hz=freq_hz, spw_warning=spw_warning,
                groups={c: [int(o) for o in oids] for c, oids in groups.items()},
                group_p95={c: float(v) for c, v in group_p95.items()}, jobs=jobs)


def run_job(plan, job, outdir):
    """One tclean (PSF + residual, niter=0) for one (field, config); returns the result row."""
    vis, datacolumn = plan["vis"], plan["datacolumn"]
    field, config, imsize = job["field"], job["config"], job["imsize"]
    prefix = os.path.join(outdir, f"test_{field}_{config}")
    row = dict(job, status="ok", psf_peak=float("nan"), resid_rms=float("nan"))
    try:
        for ext in (".psf", ".residual"):
            shutil.rmtree(prefix + ext, ignore_errors=True)
        tclean(vis=vis, imagename=prefix, field=field, observation=job["obsids"],
               datacolumn=datacolumn, specmode="cube",
               start=job["chan"], width=1, nchan=1, imsize=[imsize, imsize],
               cell=f"{job['cell']:.4f}arcsec",
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
        # Keep only the PSF and residual. The trailing '.' limits this to this job's own
        # files, so it is safe while other shards write into the same outdir.
        for fn in os.listdir(outdir):
            full = os.path.join(outdir, fn)
            if fn.startswith(os.path.basename(prefix) + ".") and \
                    not fn.endswith((".psf", ".residual")):
                shutil.rmtree(full) if os.path.isdir(full) else os.remove(full)
    return row


def write_summary(plan, rows, outdir):
    summary = os.path.join(outdir, "test_imaging_summary.txt")
    with open(summary, "w") as fh:
        fh.write(f"Test imaging for {plan['vis']}\n")
        fh.write(f"SPW {plan['spw_label']} channel {plan['chan']}/{plan['nchan']} "
                 f"({plan['freq_hz'] / 1e6:.3f} MHz), natural weighting, "
                 f"datacolumn={plan['datacolumn']}\n")
        if plan["spw_warning"]:
            fh.write(plan["spw_warning"] + "\n")
        fh.write("\n")
        fh.write("Config  obsids  p95_baseline_m\n")
        for c, o in sorted(plan["groups"].items()):
            fh.write(f"{c:>6}  {','.join(map(str, o)):>6}  {plan['group_p95'][c]:.0f}\n")
        fh.write("\n")
        fh.write("{:<16}{:>4}{:>12}{:>7}{:>10}{:>12}  status\n".format(
            "field", "cfg", "cell_arcsec", "imsize", "psf_peak", "resid_rms"))
        for r in rows:
            fh.write(f"{r['field']:<16}{r['config']:>4}{r['cell']:>12.3f}{r['imsize']:>7}"
                     f"{r['psf_peak']:>10.3f}{r['resid_rms']:>12.4g}  {r['status']}\n")
    print(f"Wrote {summary}")


def load_plan(outdir):
    with open(os.path.join(outdir, PLAN_NAME)) as fh:
        return json.load(fh)


def main():
    args = script_args("run_test_imaging.py")
    vis = args[0].rstrip("/")
    outdir = args[1] if len(args) > 1 else os.path.dirname(os.path.abspath(vis)) + "/test_images"
    # Modes are positional words, not --flags: only positional arguments are known to pass
    # through `casa -c script ...` untouched.
    mode = args[2] if len(args) > 2 else None
    os.makedirs(outdir, exist_ok=True)

    if mode is None:                       # original behaviour: everything, serially
        plan = build_plan(vis)
        rows = [run_job(plan, j, outdir) for j in plan["jobs"]]
        write_summary(plan, rows, outdir)

    elif mode == "plan":
        plan = build_plan(vis)
        with open(os.path.join(outdir, PLAN_NAME), "w") as fh:
            json.dump(plan, fh)
        print(f"Planned {len(plan['jobs'])} tclean jobs: {len(plan['groups'])} configs "
              f"{sorted(plan['groups'])}")

    elif mode == "shard":
        i, n = (int(x) for x in args[3].split("/"))
        plan = load_plan(outdir)
        mine = [(k, j) for k, j in enumerate(plan["jobs"]) if k % n == i]
        print(f"Shard {i}/{n}: {len(mine)} of {len(plan['jobs'])} jobs")
        done = {}
        for count, (k, j) in enumerate(mine, 1):
            print(f"Shard {i}/{n}: job {count}/{len(mine)}: {j['field']} {j['config']}")
            done[k] = run_job(plan, j, outdir)
            with open(os.path.join(outdir, shard_name(i)), "w") as fh:   # progress survives a crash
                json.dump(done, fh)

    elif mode == "merge":
        n = int(args[3])
        plan = load_plan(outdir)
        results = {}
        for i in range(n):
            path = os.path.join(outdir, shard_name(i))
            if os.path.exists(path):
                with open(path) as fh:
                    results.update({int(k): v for k, v in json.load(fh).items()})
        rows = []
        for k, j in enumerate(plan["jobs"]):
            rows.append(results.get(k) or dict(
                j, status=f"FAILED: no result from shard {k % n} (worker crashed or timed out)",
                psf_peak=float("nan"), resid_rms=float("nan")))
        write_summary(plan, rows, outdir)
        n_bad = sum(r["status"] != "ok" for r in rows)
        print(f"Merged {len(rows)} jobs from {n} shards; {n_bad} not ok")
        for fn in [PLAN_NAME] + [shard_name(i) for i in range(n)]:
            try:
                os.remove(os.path.join(outdir, fn))
            except OSError:
                pass
    else:
        raise SystemExit(f"Unknown mode {mode}; use plan, shard I/N or merge N")


main()
