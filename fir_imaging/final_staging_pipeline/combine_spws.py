"""Combine near-identical SPWs of a staged MS into one SPW, in place.

Some staged MSs have >1 SPW that differ only by a sub-channel frequency offset:
  - OH lines (ic1613, wlm, ngc6822 ...): a few scans on a B0D0 baseband and the
    rest on A0C0, same nchan and channel width, offset by ~0.25 channel (~2 kHz).
  - M31 hilores: 3 SPWs of the same baseband differing by ~1 Hz in width.

All SPWs used by the data are regridded onto the channel grid of the SPW that
holds the most rows (the reference) and merged with mstransform. Spectral
interpolation is 'fftshift' (fixed in CASA >= 6.7.5; run with 6.7.6).

Run with:
    casa --nologger --nogui --log2term -c combine_spws.py <ms>

After the merge, the outer edge channels are checked: regridding can leave the first
and/or last channel(s) fully flagged (no valid data in any row). If so, they are
removed with split so the staged MS has no dead edge channels (test imaging uses
channel 0, which must hold data).

MSs with a single SPW in use are left unchanged. If the SPWs differ in nchan,
differ in channel width by >1%, or are offset by >= 1 channel, the script
stops without changing the MS (these are not the small-offset case).
"""

import os
import shutil
import sys

import numpy as np
from casatasks import mstransform, split
from casatools import msmetadata, table

max_width_frac = 0.01   # allowed fractional channel-width difference
max_shift_chan = 1.0    # allowed offset of channel 0 from reference, in channels
flag_chunk_rows = 200000  # rows read at a time when counting flags per channel


def flagged_rows_per_channel(ms):
    """Return (rows in which each channel is flagged, total rows) for `ms`.

    A channel counts as flagged in a row if all correlations are flagged.
    """
    t = table()
    t.open(ms)
    nrow = t.nrows()
    counts = None
    for start in range(0, nrow, flag_chunk_rows):
        n = min(flag_chunk_rows, nrow - start)
        flag = t.getcol("FLAG", startrow=start, nrow=n)    # (ncorr, nchan, n)
        c = np.all(flag, axis=0).sum(axis=1)               # (nchan,)
        counts = c if counts is None else counts + c
    t.close()
    return counts, nrow


def script_args(script_name):
    """Arguments following `script_name` on a `casa -c script args` command line."""
    idx = max(i for i, a in enumerate(sys.argv) if a.endswith(script_name))
    return sys.argv[idx + 1:]


vis = script_args("combine_spws.py")[0].rstrip("/")

# --- SPWs actually used by the data, and the number of rows in each ---
tbl = table()
tbl.open(vis + "/DATA_DESCRIPTION")
dd_spw = tbl.getcol("SPECTRAL_WINDOW_ID")
tbl.close()
tbl.open(vis)
dd_rows = tbl.getcol("DATA_DESC_ID")
colnames = tbl.colnames()
nrow_old = tbl.nrows()
tbl.close()
rows_per_spw = {}
for dd, n in zip(*np.unique(dd_rows, return_counts=True)):
    s = int(dd_spw[dd])
    rows_per_spw[s] = rows_per_spw.get(s, 0) + int(n)
spw_ids = sorted(rows_per_spw)
print(f"SPWs in use (rows): {rows_per_spw}")

if len(spw_ids) < 2:
    print(f"{vis} has a single SPW in use; nothing to combine")
    sys.exit(0)

# --- Check the SPWs really are the small-offset case ---
ref = max(rows_per_spw, key=rows_per_spw.get)
msmd = msmetadata()
msmd.open(vis)
nchan_ref = int(msmd.nchan(ref))
freqs_ref = np.asarray(msmd.chanfreqs(ref))
width_ref = float(msmd.chanwidths(ref)[0])
for s in spw_ids:
    if s == ref:
        continue
    nchan = int(msmd.nchan(s))
    width = float(msmd.chanwidths(s)[0])
    shift = (float(msmd.chanfreqs(s)[0]) - freqs_ref[0]) / abs(width_ref)
    print(f"SPW {s} vs ref {ref}: nchan {nchan} vs {nchan_ref}, "
          f"width {width:.3f} vs {width_ref:.3f} Hz, shift {shift:+.3f} chan")
    if (nchan != nchan_ref or abs(width / width_ref - 1) > max_width_frac
            or abs(shift) >= max_shift_chan):
        msmd.close()
        raise RuntimeError(f"SPW {s} is not a small-offset copy of SPW {ref}; not combining")
msmd.close()

# --- Regrid onto the reference grid and merge into one SPW ---
result = vis[:-3] + "_combspw_tmp.ms" if vis.endswith(".ms") else vis + "_combspw_tmp"
if os.path.exists(result):
    shutil.rmtree(result)

# 'all' keeps whichever of DATA / CORRECTED_DATA are present
mstransform(vis=vis, outputvis=result, datacolumn="all",
            combinespws=True,
            regridms=True, mode="frequency", outframe="LSRK",
            nchan=nchan_ref, start=f"{freqs_ref[0]:.6f}Hz",
            width=f"{width_ref:.6f}Hz",
            interpolation="fftshift",
            keepflags=True)

# --- Verify the output before replacing the original ---
tbl.open(result)
nrow_new = tbl.nrows()
colnames_new = tbl.colnames()
tbl.close()
tbl.open(result + "/DATA_DESCRIPTION")
n_dd = tbl.nrows()
tbl.close()
msmd.open(result)
nspw_new = msmd.nspw()
nchan_new = int(msmd.nchan(0))
msmd.close()
if nrow_new != nrow_old:
    raise RuntimeError(f"Row count changed: {nrow_old} -> {nrow_new}; keeping original")
if nspw_new != 1 or n_dd != 1:
    raise RuntimeError(f"Expected 1 SPW after combine, got {nspw_new} SPWs / {n_dd} DDs; keeping original")
if nchan_new != nchan_ref:
    raise RuntimeError(f"nchan changed: {nchan_ref} -> {nchan_new}; keeping original")
for col in ("DATA", "WEIGHT_SPECTRUM"):
    if col in colnames and col not in colnames_new:
        raise RuntimeError(f"Column {col} lost in combine; keeping original")

# --- Drop fully flagged outer edge channels ---
counts, nrow_chk = flagged_rows_per_channel(result)
dead = counts == nrow_chk                     # flagged in every row
lo = 0
while lo < nchan_new and dead[lo]:
    lo += 1
hi = nchan_new - 1
while hi > lo and dead[hi]:
    hi -= 1
edge = np.r_[0:min(lo + 3, nchan_new), max(hi - 2, 0):nchan_new]
print("Flagged fraction of the outer channels after combine: "
      + ", ".join(f"ch{c}: {counts[c] / nrow_chk:.4f}" for c in sorted(set(edge.tolist()))))
if lo >= nchan_new:
    raise RuntimeError("All channels are fully flagged after combine; keeping original")
if lo > 0 or hi < nchan_new - 1:
    trimmed = result[:-3] + "_trim.ms" if result.endswith(".ms") else result + "_trim"
    if os.path.exists(trimmed):
        shutil.rmtree(trimmed)
    print(f"Fully flagged edge channels: removing {lo} leading and {nchan_new - 1 - hi} trailing "
          f"channel(s); keeping channels {lo}~{hi} ({hi - lo + 1} of {nchan_new})")
    split(vis=result, outputvis=trimmed, spw=f"0:{lo}~{hi}", datacolumn="all", keepflags=True)
    tbl.open(trimmed)
    nrow_trim = tbl.nrows()
    colnames_trim = tbl.colnames()
    tbl.close()
    msmd.open(trimmed)
    nspw_trim = msmd.nspw()
    nchan_trim = int(msmd.nchan(0))
    msmd.close()
    if nrow_trim != nrow_old or nspw_trim != 1 or nchan_trim != hi - lo + 1:
        raise RuntimeError(f"Trim check failed (rows {nrow_trim} vs {nrow_old}, SPWs {nspw_trim}, "
                           f"nchan {nchan_trim} vs {hi - lo + 1}); keeping original")
    for col in ("DATA", "WEIGHT_SPECTRUM", "CORRECTED_DATA"):
        if col in colnames_new and col not in colnames_trim:
            raise RuntimeError(f"Column {col} lost in trim; keeping original")
    shutil.rmtree(result)
    result = trimmed
else:
    print("No fully flagged edge channels; keeping all channels")

shutil.rmtree(vis)
shutil.move(result, vis)
print(f"combine_spws complete; replaced {vis}")
