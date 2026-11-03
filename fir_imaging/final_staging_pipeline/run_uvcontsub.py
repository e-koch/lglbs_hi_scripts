"""Continuum-subtract a staged MS in place with a 1st-order fit.

The fit uses the first 2 and last 2 channels of every SPW. The MS keeps its
name; the continuum-subtracted result replaces the original.

Run with:
    casa --nologger --nogui --log2term -c run_uvcontsub.py <ms>

Not run for MSs with 'nocontsub' in the name (exits without changes).
Uses the new-style CASA `uvcontsub` task (outputvis/fitspec).
"""

import os
import shutil
import sys

from casatasks import uvcontsub
from casatools import msmetadata, table



def script_args(script_name):
    """Arguments following `script_name` on a `casa -c script args` command line."""
    idx = max(i for i, a in enumerate(sys.argv) if a.endswith(script_name))
    return sys.argv[idx + 1:]


fitorder = 1
n_edge = 2   # channels used at each end of every SPW

vis = script_args("run_uvcontsub.py")[0].rstrip("/")

if "nocontsub" in os.path.basename(vis):
    print(f"{vis} is a nocontsub MS; skipping uvcontsub")
    sys.exit(0)

# --- SPWs actually used by the data, and their channel counts ---
tbl = table()
tbl.open(vis + "/DATA_DESCRIPTION")
spw_ids = sorted(set(int(s) for s in tbl.getcol("SPECTRAL_WINDOW_ID")))
tbl.close()
msmd = msmetadata()
msmd.open(vis)
nchan = {s: int(msmd.nchan(s)) for s in spw_ids}
msmd.close()

fit_parts = []
for s, n in nchan.items():
    if n <= 2 * n_edge:
        raise RuntimeError(f"SPW {s} has only {n} channels; cannot use {n_edge}+{n_edge} edge channels")
    fit_parts.append(f"{s}:0~{n_edge - 1};{n - n_edge}~{n - 1}")
fitspw = ",".join(fit_parts)
print(f"fitspw = {fitspw}")

tbl.open(vis)
datacolumn = "corrected" if "CORRECTED_DATA" in tbl.colnames() else "data"
tbl.close()
print(f"datacolumn = {datacolumn}")

result = vis[:-3] + "_contsub_tmp.ms" if vis.endswith(".ms") else vis + "_contsub_tmp"
if os.path.exists(result):
    shutil.rmtree(result)

uvcontsub(vis=vis, outputvis=result, fitspec=fitspw, fitorder=fitorder,
          datacolumn=datacolumn)

# --- Verify the output before replacing the original ---
tbl.open(result)
nrow_new = tbl.nrows()
tbl.close()
tbl.open(vis)
nrow_old = tbl.nrows()
tbl.close()
if nrow_new != nrow_old:
    raise RuntimeError(f"Row count changed: {nrow_old} -> {nrow_new}; keeping original")

shutil.rmtree(vis)
shutil.move(result, vis)
print(f"uvcontsub complete; replaced {vis}")
