"""Single-channel natural-weighted test images per field per VLA configuration.

For every field and every VLA configuration present in the MS, run tclean
(niter=0) on the middle channel of the first SPW, selecting only that
configuration's observation IDs. Only the PSF and residual are kept.

Configurations are identified per OBSERVATION_ID from the 95th-percentile uv
distance compared with absolute nominal A/B/C/D baselines
(pipeline_common.config_from_p95), not by rank.

Run with:
    casa --nologger --nogui --log2term -c run_test_imaging.py <ms> <outdir>

Writes <outdir>/test_<field>_<config>.{psf,residual} and
<outdir>/test_imaging_summary.txt.
"""

import os
import shutil
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipeline_common import good_imsize, identify_configs, script_args

from casatasks import tclean
from casatools import image as iatool
from casatools import msmetadata, table

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
chan = nchan // 2
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
    fh.write(f"{'field':<16}{'cfg':>4}{'cell\"':>9}{'imsize':>7}{'psf_peak':>10}{'resid_rms':>12}  status\n")
    for r in rows:
        fh.write(f"{r['field']:<16}{r['config']:>4}{r['cell']:>9.3f}{r['imsize']:>7}"
                 f"{r['psf_peak']:>10.3f}{r['resid_rms']:>12.4g}  {r['status']}\n")
print(f"Wrote {summary}")
