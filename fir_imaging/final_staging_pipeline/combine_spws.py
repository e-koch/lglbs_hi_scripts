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
and/or last channel(s) (almost) fully flagged. Rows from the SPW that is offset from the
reference grid have no valid data in the edge channel, so it is not always flagged in
*every* row (e.g. OH: ch0 is ~98% flagged, only the few scans on the reference SPW keep
data). Edge channels flagged in >= edge_flag_frac of the rows are removed with split so
the staged MS has no dead edge channels.

MSs with a single SPW in use are left unchanged.

Far-offset, simultaneous SPWs (e.g. M33 OH1665/OH1667 from 17B-162): the line staging
keeps every SPW that overlaps the requested line window in frequency. In 17B-162 the
OH1665 and OH1667 SPWs are 2 MHz wide and only ~2 MHz apart, so each line's window also
catches the neighbouring line's SPW. These extra SPWs are recorded in the *same scans* as
the reference SPW, are offset by many channels, and carry no data wanted for this line.
They cannot be merged (combining would merge rows and average different frequency windows),
so they are dropped with split, keeping only the reference SPW (+ any small-offset SPWs).
This is only done if every offending SPW is observed in the same observations as the
reference SPW and >= min_shared_scan_frac of its scans are shared with the reference;
otherwise the script stops without changing the MS.

If the SPWs differ in nchan, differ in channel width by >1%, or are offset by >= 1
channel and are NOT simultaneous with the reference SPW, the script stops without
changing the MS (these are not a case handled here).
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
edge_flag_frac = 0.9    # an edge channel flagged in >= this fraction of rows is dropped
min_shared_scan_frac = 0.9  # an offset SPW is "simultaneous" if >= this fraction of its
                            # (observation, scan) pairs also contain the reference SPW
scan_chunk_rows = 2000000   # rows read at a time when listing (obs, scan) per SPW


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


def spws_in_use(ms):
    """Return ({spw: number of rows}, column names, total rows) for `ms`."""
    t = table()
    t.open(ms + "/DATA_DESCRIPTION")
    dd_spw = t.getcol("SPECTRAL_WINDOW_ID")
    t.close()
    t.open(ms)
    dd_rows = t.getcol("DATA_DESC_ID")
    cols = t.colnames()
    nrow = t.nrows()
    t.close()
    rows = {}
    for dd, n in zip(*np.unique(dd_rows, return_counts=True)):
        s = int(dd_spw[dd])
        rows[s] = rows.get(s, 0) + int(n)
    return rows, cols, nrow


def not_small_offset(ms, spw_ids, ref):
    """Return the SPWs in `spw_ids` that are not small-offset copies of `ref`."""
    md = msmetadata()
    md.open(ms)
    nchan_ref = int(md.nchan(ref))
    freq0_ref = float(md.chanfreqs(ref)[0])
    width_ref = float(md.chanwidths(ref)[0])
    bad = []
    for s in spw_ids:
        if s == ref:
            continue
        nchan = int(md.nchan(s))
        width = float(md.chanwidths(s)[0])
        shift = (float(md.chanfreqs(s)[0]) - freq0_ref) / abs(width_ref)
        print(f"SPW {s} vs ref {ref}: nchan {nchan} vs {nchan_ref}, "
              f"width {width:.3f} vs {width_ref:.3f} Hz, shift {shift:+.3f} chan")
        if (nchan != nchan_ref or abs(width / width_ref - 1) > max_width_frac
                or abs(shift) >= max_shift_chan):
            bad.append(s)
    md.close()
    return bad


def scans_per_spw(ms):
    """Return {spw: set of (observation_id, scan_number)} for the rows of `ms`."""
    t = table()
    t.open(ms + "/DATA_DESCRIPTION")
    dd_spw = t.getcol("SPECTRAL_WINDOW_ID")
    t.close()
    t.open(ms)
    nrow = t.nrows()
    pairs = {}
    for start in range(0, nrow, scan_chunk_rows):
        n = min(scan_chunk_rows, nrow - start)
        dd = t.getcol("DATA_DESC_ID", startrow=start, nrow=n)
        obs = t.getcol("OBSERVATION_ID", startrow=start, nrow=n)
        scan = t.getcol("SCAN_NUMBER", startrow=start, nrow=n)
        uniq = np.unique(np.stack([dd, obs, scan]), axis=1)
        for d, o, sc in uniq.T:
            pairs.setdefault(int(dd_spw[int(d)]), set()).add((int(o), int(sc)))
    t.close()
    return pairs


def simultaneous_with_ref(spw_scans, s, ref):
    """True if SPW `s` is only seen in observations that also contain `ref`, and
    >= min_shared_scan_frac of its (observation, scan) pairs also contain `ref`.
    Returns (ok, n_scans_of_s, n_shared_scans, observations_of_s_without_ref)."""
    s_scans = spw_scans[s]
    shared = s_scans & spw_scans[ref]
    obs_without_ref = {o for o, _ in s_scans} - {o for o, _ in spw_scans[ref]}
    frac = len(shared) / max(len(s_scans), 1)
    return (not obs_without_ref and frac >= min_shared_scan_frac,
            len(s_scans), len(shared), sorted(obs_without_ref))


def drop_spws_keep(ms, keep_spws):
    """Replace `ms` (in place) by a copy that contains only `keep_spws` (split)."""
    kept = f"{ms[:-3]}_keepspw_tmp.ms" if ms.endswith(".ms") else ms + "_keepspw_tmp"
    if os.path.exists(kept):
        shutil.rmtree(kept)
    rows_old, cols_old, nrow_all = spws_in_use(ms)
    nrow_expect = sum(rows_old[s] for s in keep_spws)
    split(vis=ms, outputvis=kept, spw=",".join(str(s) for s in keep_spws),
          datacolumn="all", keepflags=True)
    rows_new, cols_new, nrow_new = spws_in_use(kept)
    if nrow_new != nrow_expect or len(rows_new) != len(keep_spws):
        raise RuntimeError(f"SPW selection check failed (rows {nrow_new} vs {nrow_expect}, "
                           f"SPWs {len(rows_new)} vs {len(keep_spws)}); keeping original")
    for col in ("DATA", "WEIGHT_SPECTRUM", "CORRECTED_DATA"):
        if col in cols_old and col not in cols_new:
            raise RuntimeError(f"Column {col} lost in SPW selection; keeping original")
    shutil.rmtree(ms)
    shutil.move(kept, ms)
    print(f"Dropped {nrow_all - nrow_new} rows ({len(rows_old) - len(rows_new)} SPWs); "
          f"replaced {ms}")


vis = script_args("combine_spws.py")[0].rstrip("/")

# --- SPWs actually used by the data, and the number of rows in each ---
rows_per_spw, colnames, nrow_old = spws_in_use(vis)
spw_ids = sorted(rows_per_spw)
print(f"SPWs in use (rows): {rows_per_spw}")

if len(spw_ids) < 2:
    print(f"{vis} has a single SPW in use; nothing to combine")
    sys.exit(0)

# --- Check the SPWs really are the small-offset case ---
ref = max(rows_per_spw, key=rows_per_spw.get)
bad = not_small_offset(vis, spw_ids, ref)

if bad:
    # Not small-offset copies. Only OK if they are simultaneous with the reference
    # (extra windows recorded in the same scans): drop them instead of merging.
    spw_scans = scans_per_spw(vis)
    not_sim = []
    for s in bad:
        ok, n_s, n_shared, obs_missing = simultaneous_with_ref(spw_scans, s, ref)
        print(f"SPW {s}: {n_s} (obs, scan) pairs, {n_shared} shared with ref SPW {ref}"
              + (f"; observations without ref: {obs_missing}" if obs_missing else ""))
        if not ok:
            not_sim.append(s)
    if not_sim:
        raise RuntimeError(f"SPW(s) {not_sim} are not small-offset copies of SPW {ref} and are "
                           f"not simultaneous with it; not combining")
    print(f"Dropping far-offset SPWs {bad} (recorded in the same scans as SPW {ref}); "
          f"keeping SPWs {[s for s in spw_ids if s not in bad]}")
    drop_spws_keep(vis, [s for s in spw_ids if s not in bad])

    # Re-read the (renumbered) SPWs in the reduced MS
    rows_per_spw, colnames, nrow_old = spws_in_use(vis)
    spw_ids = sorted(rows_per_spw)
    print(f"SPWs in use after dropping (rows): {rows_per_spw}")
    if len(spw_ids) < 2:
        print(f"{vis} has a single SPW in use; nothing to combine")
        sys.exit(0)
    ref = max(rows_per_spw, key=rows_per_spw.get)
    if not_small_offset(vis, spw_ids, ref):
        raise RuntimeError("Remaining SPWs are still not small-offset copies; not combining")

msmd = msmetadata()
msmd.open(vis)
nchan_ref = int(msmd.nchan(ref))
freqs_ref = np.asarray(msmd.chanfreqs(ref))
width_ref = float(msmd.chanwidths(ref)[0])
msmd.close()
tbl = table()

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

# --- Drop (almost) fully flagged outer edge channels ---
counts, nrow_chk = flagged_rows_per_channel(result)
dead = counts >= edge_flag_frac * nrow_chk    # flagged in (nearly) every row
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
    print(f"Edge channels flagged in >= {edge_flag_frac:.0%} of rows: removing {lo} leading and {nchan_new - 1 - hi} trailing "
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
    print(f"No edge channels flagged in >= {edge_flag_frac:.0%} of rows; keeping all channels")

shutil.rmtree(vis)
shutil.move(result, vis)
print(f"combine_spws complete; replaced {vis}")
