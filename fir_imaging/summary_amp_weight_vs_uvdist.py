"""
Summary diagnostic plots of amplitude and weight vs uv-distance for a CASA MS,
plus inferred VLA array-configuration groupings.

The MS does not record which VLA configuration (A/B/C/D) each OBSERVATION_ID
("obsid") belongs to, so configs are inferred by clustering obsids on the
95th-percentile baseline length of their visibilities. Consecutive VLA configs
differ in maximum baseline length by a factor of ~3.3 (A:36 km, B:11 km,
C:3.4 km, D:1 km, i.e. ~0.5 dex steps), so a gap in log10(p95 baseline) cleanly
separates them. Any number of groups (>= 1) is supported — e.g. an MS combining
only B+C will yield two groups, not four.

All visibility data are streamed from the MS in row chunks and reduced
on-the-fly to per-(obsid, uv-distance bin) running sums (count, sum, sum of
squares), so memory use stays low and roughly constant regardless of MS size.

Run with:
    casa --nologger --nogui --log2term -c summary_amp_weight_vs_uvdist.py <vis>

Outputs (written next to the MS):
  <vis>_obsid_config_groups.txt  : obsid -> inferred config table, and
                                   per-config amplitude/weight summary table
  <vis>_amp_vs_uvdist.png        : binned amplitude vs uv-distance, one curve
                                   per inferred config group
  <vis>_weight_vs_uvdist.png     : binned weight vs uv-distance, one curve
                                   per inferred config group
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import shutil
import sys

vis = sys.argv[-1]

chunk_size        = 200_000  # rows read from the MS per chunk
target_chan       = 1        # channel sampled for amplitude/weight (keeps memory low)
n_uvdist_bins     = 40       # log-spaced uv-distance bins for the binned statistics
uvdist_min        = 20.0     # m   (below the VLA's shortest spacing)
uvdist_max        = 4.0e4    # m   (above the VLA A-config max baseline, ~36 km)
log_gap_threshold = 0.3      # dex; obsids whose log10(p95 baseline) differ by more
                             # than this are placed in different config groups
                             # (VLA configs step by ~0.5 dex; 0.3 dex leaves margin
                             # for epoch-to-epoch scatter within a config)

apply_nyquist_weights = True  # If True: copy the MS, run initweights(wtmode='nyq') on the
                               # copy, run the analysis on that re-weighted copy, then delete
                               # the copy. Lets you compare amp/weight summaries against
                               # Nyquist-based weights without altering the original MS.
                               # Output filenames get a "_nyq" suffix to keep them distinct
                               # from a plain run on the same MS.

vis_base = os.path.splitext(vis.rstrip('/').rstrip(os.sep))[0]

ms_copy     = None
weight_note = ""
if apply_nyquist_weights:
    ms_copy = f"{vis.rstrip('/').rstrip(os.sep)}_nyq_copy"
    if os.path.exists(ms_copy):
        shutil.rmtree(ms_copy)
    print(f"Copying MS to {ms_copy} ...")
    shutil.copytree(vis, ms_copy)
    print("Running initweights(wtmode='nyq') on the copy ...")
    initweights(vis=ms_copy, wtmode='nyq')
    vis_to_use  = ms_copy
    out_base    = f"{vis_base}_nyq"
    weight_note = " — weights re-initialized with initweights(wtmode='nyq')"
else:
    vis_to_use = vis
    out_base   = vis_base

# --- Inspect MS structure once ---
tb.open(vis_to_use)
colnames            = tb.colnames()
data_col            = 'CORRECTED_DATA' if 'CORRECTED_DATA' in colnames else 'DATA'
use_weight_spectrum = 'WEIGHT_SPECTRUM' in colnames
nrows               = tb.nrows()
tb.close()

print(f"MS            : {vis}")
if apply_nyquist_weights:
    print(f"Analysis copy : {ms_copy}  (initweights wtmode='nyq')")
print(f"Rows          : {nrows}")
print(f"Data column   : {data_col}")
print(f"Weight column : {'WEIGHT_SPECTRUM' if use_weight_spectrum else 'WEIGHT'}")
print(f"Channel       : {target_chan}")
print(f"Chunk size    : {chunk_size} rows")

bin_edges   = np.logspace(np.log10(uvdist_min), np.log10(uvdist_max), n_uvdist_bins + 1)
bin_centers = np.sqrt(bin_edges[:-1] * bin_edges[1:])
n_bins      = n_uvdist_bins


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _new_accumulator():
    return {
        'count':     np.zeros(n_bins, dtype=np.int64),
        'amp_sum':   np.zeros(n_bins),
        'amp_sumsq': np.zeros(n_bins),
        'wt_sum':    np.zeros(n_bins),
        'wt_sumsq':  np.zeros(n_bins),
    }


def _read_chunk(start, nrow):
    """Read one row chunk; return (obsid, uvdist, amp, weight) for unflagged rows.

    Only `target_chan` is read from DATA/WEIGHT_SPECTRUM/FLAG to keep the
    per-chunk memory footprint small and independent of the channel count.
    """
    obsid  = tb.getcol('OBSERVATION_ID', startrow=start, nrow=nrow)
    uvw    = tb.getcol('UVW', startrow=start, nrow=nrow)              # (3, nrow)
    uvdist = np.sqrt(uvw[0]**2 + uvw[1]**2)

    data_slice = tb.getcolslice(data_col, blc=[0, target_chan], trc=[-1, target_chan], incr=[],
                                startrow=start, nrow=nrow)            # (npol, 1, nrow)
    flag_slice = tb.getcolslice('FLAG', blc=[0, target_chan], trc=[-1, target_chan], incr=[],
                                startrow=start, nrow=nrow)

    amp     = np.abs(data_slice[:, 0, :]).mean(axis=0)                # (nrow,)
    flagged = flag_slice[:, 0, :].any(axis=0)

    if use_weight_spectrum:
        w_slice = tb.getcolslice('WEIGHT_SPECTRUM', blc=[0, target_chan], trc=[-1, target_chan], incr=[],
                                 startrow=start, nrow=nrow)
        weight = w_slice[:, 0, :].mean(axis=0)
    else:
        weight = tb.getcol('WEIGHT', startrow=start, nrow=nrow).mean(axis=0)

    good = (~flagged) & np.isfinite(amp) & np.isfinite(weight) & (uvdist > 0)
    return obsid[good], uvdist[good], amp[good], weight[good]


def _accumulate(acc, obsid, uvdist, amp, weight):
    """Bin by uv-distance and add to the running sums for each obsid present."""
    bin_idx = np.clip(np.digitize(uvdist, bin_edges) - 1, 0, n_bins - 1)
    for oid in np.unique(obsid):
        sel = obsid == oid
        a   = acc.setdefault(int(oid), _new_accumulator())
        idx = bin_idx[sel]
        np.add.at(a['count'],     idx, 1)
        np.add.at(a['amp_sum'],   idx, amp[sel])
        np.add.at(a['amp_sumsq'], idx, amp[sel]**2)
        np.add.at(a['wt_sum'],    idx, weight[sel])
        np.add.at(a['wt_sumsq'],  idx, weight[sel]**2)


def _percentile_from_hist(counts, pct):
    """Linearly-interpolated percentile (0-100) from a binned histogram of uv-distance."""
    total = counts.sum()
    if total == 0:
        return np.nan
    cum    = np.cumsum(counts) / total
    target = pct / 100.0
    i      = min(int(np.searchsorted(cum, target)), n_bins - 1)
    cum_lo = cum[i - 1] if i > 0 else 0.0
    cum_hi = cum[i]
    frac   = 0.0 if cum_hi == cum_lo else (target - cum_lo) / (cum_hi - cum_lo)
    return bin_edges[i] + frac * (bin_edges[i + 1] - bin_edges[i])


def _cluster_by_log_gap(value_by_key, gap):
    """Sort keys by value (desc); start a new group wherever the log10 gap exceeds `gap`."""
    order  = sorted(value_by_key, key=value_by_key.get, reverse=True)
    groups = [[order[0]]]
    for prev, cur in zip(order, order[1:]):
        if np.log10(value_by_key[prev]) - np.log10(value_by_key[cur]) <= gap:
            groups[-1].append(cur)
        else:
            groups.append([cur])
    return groups


def _config_label(i):
    letters = 'ABCD'
    return letters[i] if i < len(letters) else f'GRP{i + 1}'


def _bin_mean_std(sum_, sumsq_, count_):
    with np.errstate(invalid='ignore', divide='ignore'):
        mean = np.where(count_ > 0, sum_ / np.maximum(count_, 1), np.nan)
        var  = np.where(count_ > 0, sumsq_ / np.maximum(count_, 1) - mean**2, np.nan)
    return mean, np.sqrt(np.clip(var, 0, None))


def _print_table(title, headers, rows):
    print(f"\n--- {title} ---")
    if not rows:
        print("(no rows)")
        return
    widths = [max(len(h), *(len(str(r[c])) for r in rows)) for c, h in enumerate(headers)]
    line = '  '.join(h.rjust(w) for h, w in zip(headers, widths))
    print(line)
    print('-' * len(line))
    for r in rows:
        print('  '.join(str(v).rjust(w) for v, w in zip(r, widths)))


def _write_table(fh, title, headers, rows):
    fh.write(f"{title}\n")
    if not rows:
        fh.write("(no rows)\n\n")
        return
    widths = [max(len(h), *(len(str(r[c])) for r in rows)) for c, h in enumerate(headers)]
    line = '  '.join(h.rjust(w) for h, w in zip(headers, widths))
    fh.write(line + '\n')
    fh.write('-' * len(line) + '\n')
    for r in rows:
        fh.write('  '.join(str(v).rjust(w) for v, w in zip(r, widths)) + '\n')
    fh.write('\n')


# ---------------------------------------------------------------------------
# Stream the MS in chunks, accumulating per-obsid binned statistics
# ---------------------------------------------------------------------------

acc = {}   # obsid -> accumulator dict (count, amp_sum, amp_sumsq, wt_sum, wt_sumsq)
tb.open(vis_to_use)
n_chunks = int(np.ceil(nrows / chunk_size))
for ichunk in range(n_chunks):
    start = ichunk * chunk_size
    nrow  = min(chunk_size, nrows - start)
    obsid, uvdist, amp, weight = _read_chunk(start, nrow)
    _accumulate(acc, obsid, uvdist, amp, weight)
    print(f"  Chunk {ichunk + 1}/{n_chunks}: rows {start}-{start + nrow - 1} "
          f"({len(obsid)} unflagged)")
tb.close()

obsids = sorted(acc)
print(f"\nFound {len(obsids)} observation ID(s): {obsids}")


# ---------------------------------------------------------------------------
# Cluster obsids into inferred VLA configurations by 95th-percentile baseline
# ---------------------------------------------------------------------------

p95_baseline = {oid: _percentile_from_hist(acc[oid]['count'], 95.0) for oid in obsids}

usable_obsids = [oid for oid in obsids if np.isfinite(p95_baseline[oid])]
skipped       = sorted(set(obsids) - set(usable_obsids))
if skipped:
    print(f"Warning: no usable (unflagged) rows for obsid(s) {skipped}; excluding from grouping")

groups = _cluster_by_log_gap({oid: p95_baseline[oid] for oid in usable_obsids}, log_gap_threshold)

group_obsids = {_config_label(i): sorted(grp) for i, grp in enumerate(groups)}

print(f"\nInferred {len(groups)} configuration group(s) from log10(p95 baseline) "
      f"clustering (gap threshold = {log_gap_threshold} dex):")
for label, grp in group_obsids.items():
    p95s = [p95_baseline[oid] for oid in grp]
    print(f"  Config {label}: obsid(s) {grp}, p95 baseline ~ {min(p95s):.0f}-{max(p95s):.0f} m")

# --- Table 1: obsid -> inferred configuration ---
rows1 = [(oid, label, f"{p95_baseline[oid]:.0f}", int(acc[oid]['count'].sum()))
         for label, grp in group_obsids.items() for oid in grp]
headers1 = ('obsid', 'config', 'p95_baseline_m', 'n_rows')
_print_table("Obsid -> inferred VLA configuration", headers1, rows1)


# ---------------------------------------------------------------------------
# Aggregate per-config statistics by summing the per-obsid binned accumulators
# ---------------------------------------------------------------------------

group_acc = {}
for label, grp in group_obsids.items():
    g = _new_accumulator()
    for oid in grp:
        for key in g:
            g[key] = g[key] + acc[oid][key]
    group_acc[label] = g

# --- Table 2: per-config amplitude & weight summary ---
rows2 = []
for label, g in group_acc.items():
    n = int(g['count'].sum())
    amp_m, amp_s = _bin_mean_std(np.array([g['amp_sum'].sum()]), np.array([g['amp_sumsq'].sum()]), np.array([n]))
    wt_m, wt_s   = _bin_mean_std(np.array([g['wt_sum'].sum()]),  np.array([g['wt_sumsq'].sum()]),  np.array([n]))
    rows2.append((label, len(group_obsids[label]), n,
                  f"{amp_m[0]:.4g} +/- {amp_s[0]:.2g}",
                  f"{wt_m[0]:.4g} +/- {wt_s[0]:.2g}"))
headers2 = ('config', 'n_obsids', 'n_rows', 'amplitude (mean+/-std)', 'weight (mean+/-std)')
_print_table("Per-configuration amplitude & weight summary", headers2, rows2)

# --- Save both tables to a text file ---
txt_path = f"{out_base}_obsid_config_groups.txt"
with open(txt_path, 'w') as fh:
    fh.write(f"# {vis}\n")
    if apply_nyquist_weights:
        fh.write(f"# Weights re-initialized via initweights(wtmode='nyq') on a temporary "
                 f"copy ({ms_copy})\n")
    fh.write(f"# Channel sampled       : {target_chan}\n")
    fh.write(f"# Config groups inferred from clustering log10(95th-pct baseline length)\n")
    fh.write(f"# log10 gap threshold   : {log_gap_threshold} dex\n\n")
    _write_table(fh, "Obsid -> inferred VLA configuration", headers1, rows1)
    _write_table(fh, "Per-configuration amplitude & weight summary", headers2, rows2)
print(f"\nSaved tables: {txt_path}")


# ---------------------------------------------------------------------------
# Summary plots: binned amplitude / weight vs uv-distance, per config group
# ---------------------------------------------------------------------------

def _plot_vs_uvdist(sum_key, sumsq_key, ylabel, log_y, out_path):
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = plt.cm.viridis(np.linspace(0, 1, max(len(group_acc), 2)))
    for color, (label, g) in zip(colors, group_acc.items()):
        mean, std = _bin_mean_std(g[sum_key], g[sumsq_key], g['count'])
        valid = g['count'] > 0
        if not valid.any():
            continue
        ax.plot(bin_centers[valid], mean[valid], '-o', ms=4, color=color,
                label=f"Config {label} ({len(group_obsids[label])} obsid(s), "
                      f"{int(g['count'].sum())} rows)")
        ax.fill_between(bin_centers[valid], (mean - std)[valid], (mean + std)[valid],
                        color=color, alpha=0.2)
    ax.set_xscale('log')
    if log_y:
        ax.set_yscale('log')
    ax.set_xlabel('uv-distance [m]')
    ax.set_ylabel(ylabel)
    ax.set_title(f'{vis_base}{weight_note}\n'
                 f'channel {target_chan} — binned mean +/- std per inferred config')
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved plot: {out_path}")


_plot_vs_uvdist('amp_sum', 'amp_sumsq', 'Amplitude', False, f"{out_base}_amp_vs_uvdist.png")
_plot_vs_uvdist('wt_sum',  'wt_sumsq',  'Weight',    True,  f"{out_base}_weight_vs_uvdist.png")

# --- Clean up the temporary re-weighted MS copy ---
if apply_nyquist_weights and ms_copy is not None:
    print(f"\nRemoving temporary MS copy {ms_copy} ...")
    shutil.rmtree(ms_copy)
