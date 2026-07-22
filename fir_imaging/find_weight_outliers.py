"""
Find severe outliers in the visibility weight distribution of a CASA MS.

Loads a single channel from each SPW to keep memory usage low.
Prints a ranked list of outlier baselines with scan and antenna labels
suitable for manual flagging with flagdata(). Per SPW, saves:
  - <vis>_spw<N>_weight_outliers.npy  : scan/antenna dicts
  - <vis>_spw<N>_flagcmds.txt         : flagdata mode='list' commands
  - <vis>_spw<N>_<ant>_scans<X>_weights.png : weight vs scan whisker plots

Run with:
    casa --nologger --nogui --log2term -c find_weight_outliers.py <vis>

Load saved output:
    d = np.load('your_spw0_weight_outliers.npy', allow_pickle=True).item()
    # d['inf_weights']     -> {scan: {'all': [...], 'hub': [...], 'hub_fractions': {...}}}
    # d['finite_outliers'] -> {scan: {'all': [...], 'hub': [...], 'hub_fractions': {...}}}
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from datetime import datetime, timedelta

_MJD_EPOCH = datetime(1858, 11, 17)

vis             = sys.argv[-1]
target_chan     = 1      # which channel to sample for weight statistics
threshold_sigma = 100.0  # flag rows where weight > median + N * scaled-MAD
hub_threshold   = 0.9    # fraction of scan antennas a hub must have outlier baselines with
scan_context    = 10     # number of scans either side of an outlier group to plot

_ = listobs(vis, listfile=f"{vis}.listobs.txt", overwrite=True)

# --- Get observation ID and date ---
tb.open(vis + '/OBSERVATION')
obs_project    = tb.getcol('PROJECT')
obs_time_range = tb.getcol('TIME_RANGE')   # shape (2, nobs): [start, end] MJD seconds
tb.close()
obs_id   = str(obs_project[0]) if len(obs_project) else 'unknown'
obs_date = (_MJD_EPOCH + timedelta(seconds=float(obs_time_range[0, 0]))).strftime('%d-%m-%Y')
print(f"Observation ID : {obs_id}")
print(f"Observation date: {obs_date}")

# --- Get antenna names ---
tb.open(vis + '/ANTENNA')
ant_names = tb.getcol('NAME')
tb.close()

# --- Get DATA_DESC_ID -> SPW mapping ---
tb.open(vis + '/DATA_DESCRIPTION')
spw_ids = tb.getcol('SPECTRAL_WINDOW_ID')
tb.close()
n_spws = len(spw_ids)
print(f"DATA_DESC_ID -> SPW mapping: {list(enumerate(spw_ids))}")

# --- Check once whether WEIGHT_SPECTRUM exists ---
tb.open(vis)
use_weight_spectrum = 'WEIGHT_SPECTRUM' in tb.colnames()
tb.close()
print(f"Using {'WEIGHT_SPECTRUM' if use_weight_spectrum else 'WEIGHT'}")

vis_base = os.path.splitext(vis.rstrip('/').rstrip(os.sep))[0]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _scan_antenna_dict(mask, scans, ant1, ant2):
    """Return {scan: {'all': [...], 'hub': [...], 'hub_fractions': {...}}}.

    'all'           — every antenna on any outlier baseline in this scan.
    'hub'           — antennas whose outlier baselines cover >= hub_threshold
                      fraction of all other antennas present in the scan.
                      These are the likely root-cause antennas.
    'hub_fractions' — {antenna: fraction} for every antenna in 'all', for
                      inspection when tuning hub_threshold.
    """
    # Total antennas present per scan (from all rows, not just outliers)
    scan_all_ants = {}
    for i in range(len(scans)):
        scan = int(scans[i])
        if scan not in scan_all_ants:
            scan_all_ants[scan] = set()
        scan_all_ants[scan].add(ant_names[ant1[i]])
        scan_all_ants[scan].add(ant_names[ant2[i]])

    # Collect outlier baselines per scan
    outlier_baselines = {}
    for i in np.where(mask)[0]:
        scan = int(scans[i])
        if scan not in outlier_baselines:
            outlier_baselines[scan] = []
        outlier_baselines[scan].append((ant_names[ant1[i]], ant_names[ant2[i]]))

    result = {}
    for scan, bl_list in sorted(outlier_baselines.items()):
        all_ants = sorted(set(a for pair in bl_list for a in pair))
        n_total  = len(scan_all_ants[scan])

        # Count distinct outlier partners per antenna
        partners = {}
        for a, b in bl_list:
            partners.setdefault(a, set()).add(b)
            partners.setdefault(b, set()).add(a)

        # Fraction = distinct outlier partners / (total ants in scan - 1)
        fractions = {ant: len(p) / (n_total - 1) for ant, p in partners.items()}
        hub_ants  = sorted(ant for ant, frac in fractions.items()
                           if frac >= hub_threshold)

        result[scan] = {'all': all_ants, 'hub': hub_ants, 'hub_fractions': fractions}
    return result


def _print_block(label, mask, w, scans, ant1, ant2, hub_ants):
    """Print outlier rows restricted to those involving a hub antenna."""
    hub_set = set(hub_ants)
    hub_row = np.array([(ant_names[ant1[i]] in hub_set or ant_names[ant2[i]] in hub_set)
                        for i in range(len(scans))])
    idx = np.where(mask & hub_row)[0]
    if len(idx) == 0:
        return
    idx = idx[np.argsort(w[idx])[::-1]]
    print()
    print(f"--- {label} ---")
    print(f"{'Rank':>4}  {'Weight':>10}  {'Scan':>6}  {'Ant1':>6}  {'Ant2':>6}  {'Baseline':>12}")
    print("-" * 60)
    for rank, i in enumerate(idx, 1):
        a1 = ant_names[ant1[i]]
        a2 = ant_names[ant2[i]]
        weight_str = '       inf' if np.isinf(w[i]) else f"{w[i]:>10.2f}"
        print(f"{rank:>4}  {weight_str}  {scans[i]:>6}  {a1:>6}  {a2:>6}  {a1}-{a2}")
    bad_scans = np.unique(scans[idx])
    print()
    print("Scans    :", sorted(bad_scans.tolist()))
    print("Hub ants :", sorted(hub_set))
    print()
    scan_str = ','.join(str(s) for s in sorted(bad_scans.tolist()))
    ant_str  = ','.join(sorted(hub_set))
    print(f"  flagdata(vis='{vis}', mode='manual', scan='{scan_str}', antenna='{ant_str}')")


def _group_sequential_scans(scan_list, all_scans_sorted):
    """Group scans that are adjacent in all_scans_sorted into contiguous runs.

    Returns a list of lists, each inner list being one sequential group.
    """
    if not scan_list:
        return []
    positions = sorted(int(np.searchsorted(all_scans_sorted, s)) for s in scan_list)
    groups = [[positions[0]]]
    for pos in positions[1:]:
        if pos == groups[-1][-1] + 1:
            groups[-1].append(pos)
        else:
            groups.append([pos])
    return [[int(all_scans_sorted[p]) for p in grp] for grp in groups]


def _plot_hub_weights(vis_base, spw, hub_ant, scan_group,
                      all_scans_sorted, ant1_names, ant2_names, scans, w):
    """Whisker plot of weights vs scan for ±scan_context scans around scan_group.

    Box/whisker shows all non-hub-antenna baselines per scan.
    Hub antenna baselines are overlaid as individual scatter points.
    Infinite weights for the hub antenna are shown as separate markers at
    the top of the plot.
    """
    # ±scan_context window around the group extent
    pos_lo = int(np.searchsorted(all_scans_sorted, scan_group[0]))
    pos_hi = int(np.searchsorted(all_scans_sorted, scan_group[-1]))
    lo = max(0, pos_lo - scan_context)
    hi = min(len(all_scans_sorted), pos_hi + scan_context + 1)
    plot_scans = all_scans_sorted[lo:hi]

    is_hub = (ant1_names == hub_ant) | (ant2_names == hub_ant)

    box_data  = []   # per scan: finite non-hub weights
    hub_x     = []   # scatter x position (1-based index)
    hub_y     = []   # scatter y (finite hub weights)
    inf_x     = []   # x positions of inf hub weights

    for pos, scan in enumerate(plot_scans, start=1):
        scan_mask = scans == scan
        other_w = w[scan_mask & ~is_hub]
        hub_w   = w[scan_mask &  is_hub]

        box_data.append(other_w[np.isfinite(other_w)] if np.isfinite(other_w).any()
                        else np.array([np.nan]))

        finite_hub = hub_w[np.isfinite(hub_w)]
        hub_x.extend([pos] * len(finite_hub))
        hub_y.extend(finite_hub.tolist())

        if np.isinf(hub_w).any():
            inf_x.append(pos)

    fig, ax = plt.subplots(figsize=(max(8, len(plot_scans) * 0.45), 5))

    ax.boxplot(box_data,
               positions=np.arange(1, len(plot_scans) + 1),
               widths=0.6, sym='', patch_artist=True,
               boxprops=dict(facecolor='steelblue', alpha=0.5),
               medianprops=dict(color='navy'))

    if hub_x:
        ax.scatter(hub_x, hub_y, color='red', marker='*', s=100, zorder=5,
                   label=f'{hub_ant} (finite)')

    if inf_x:
        # Draw inf markers just above the plot top; will be clipped to axes
        ylim_top = ax.get_ylim()[1]
        ax.scatter(inf_x, [ylim_top] * len(inf_x),
                   color='darkred', marker='^', s=120, zorder=5, clip_on=False,
                   label=f'{hub_ant} (inf)')

    # Shade outlier scans
    for scan in scan_group:
        matches = np.where(plot_scans == scan)[0]
        if len(matches):
            pos = matches[0] + 1
            ax.axvspan(pos - 0.5, pos + 0.5, color='red', alpha=0.15, zorder=0)

    ax.set_xticks(np.arange(1, len(plot_scans) + 1))
    ax.set_xticklabels([str(s) for s in plot_scans],
                       rotation=45, ha='right', fontsize=8)
    ax.set_xlabel('Scan')
    ax.set_ylabel('Weight')
    ax.set_title(f'SPW {spw} — hub antenna: {hub_ant}\n'
                 f'Outlier scans: {scan_group}  (shaded)')
    ax.legend(loc='upper left')
    fig.tight_layout()

    scan_label = (f"{scan_group[0]}-{scan_group[-1]}"
                  if len(scan_group) > 1 else str(scan_group[0]))
    out_path = f"{vis_base}_spw{spw}_{hub_ant}_scans{scan_label}_weights.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  Saved plot: {out_path}")


# ---------------------------------------------------------------------------
# Main loop over SPWs
# ---------------------------------------------------------------------------

for ddid in range(n_spws):
    spw = spw_ids[ddid]
    print(f"\n{'='*60}")
    print(f"DATA_DESC_ID {ddid}  (SPW {spw})")
    print(f"{'='*60}")

    tb.open(vis)
    subt = tb.query(f'DATA_DESC_ID == {ddid}')

    ant1  = subt.getcol('ANTENNA1')
    ant2  = subt.getcol('ANTENNA2')
    scans = subt.getcol('SCAN_NUMBER')

    if use_weight_spectrum:
        ws = subt.getcolslice('WEIGHT_SPECTRUM',
                               blc=[0, target_chan],
                               trc=[-1, target_chan])   # shape (npol, 1, nrow)
        w = ws[:, 0, :].max(axis=0)                     # shape (nrow,)
    else:
        w = subt.getcol('WEIGHT').max(axis=0)            # shape (nrow,)

    subt.close()
    tb.close()

    assert len(ant1) == len(ant2) == len(scans) == len(w)
    print(f"Loaded {len(w)} rows")

    # Pre-compute antenna name arrays once per SPW (used in plotting)
    ant1_names = np.array([ant_names[a] for a in ant1])
    ant2_names = np.array([ant_names[a] for a in ant2])
    all_scans_sorted = np.array(sorted(np.unique(scans).tolist()))

    # --- Separate infs from finite weights ---
    inf_mask    = np.isinf(w)
    finite_mask = ~inf_mask
    w_finite    = w[finite_mask]

    # --- Robust outlier detection on finite values only ---
    med = np.median(w_finite)
    mad = np.median(np.abs(w_finite - med))
    robust_std = 1.4826 * mad
    upper = med + threshold_sigma * robust_std

    outlier_mask = finite_mask & (w > upper)
    n_inf      = inf_mask.sum()
    n_outliers = outlier_mask.sum()

    print(f"Median weight      : {med:.4f}")
    print(f"MAD (scaled)       : {robust_std:.4f}")
    print(f"Threshold          : {upper:.4f}  ({threshold_sigma}-sigma)")
    print(f"Inf rows           : {n_inf} / {len(w)}")
    print(f"Finite outlier rows: {n_outliers} / {len(w)}")

    # --- Build scan/antenna dicts (needed before printing to get hub ants) ---
    output = {
        'obs_id':          obs_id,
        'obs_date':        obs_date,
        'inf_weights':     _scan_antenna_dict(inf_mask,     scans, ant1, ant2),
        'finite_outliers': _scan_antenna_dict(outlier_mask, scans, ant1, ant2),
    }

    if n_inf > 0:
        inf_hub_ants = sorted({ant for info in output['inf_weights'].values()
                                for ant in info['hub']})
        print(f"Found {n_inf} infinite weights")
        if inf_hub_ants:
            _print_block("Infinite weights — hub rows", inf_mask, w, scans, ant1, ant2,
                         inf_hub_ants)

    if n_outliers > 0:
        outlier_hub_ants = sorted({ant for info in output['finite_outliers'].values()
                                    for ant in info['hub']})
        if outlier_hub_ants:
            _print_block(f"Finite outliers (>{threshold_sigma}-sigma) — hub rows",
                         outlier_mask, w, scans, ant1, ant2, outlier_hub_ants)
        else:
            print("No hub antennas identified among finite outliers.")

    # --- Save .npy ---
    npy_path = f"{vis_base}_spw{spw}_weight_outliers.npy"
    np.save(npy_path, output)
    print(f"\nSaved dict: {npy_path}")

    # --- Build flagdata mode='list' commands from hub antennas ---
    flag_lines = []
    for category, label in [('inf_weights', 'inf'), ('finite_outliers', 'outlier')]:
        for scan, info in output[category].items():
            for hub_ant in info['hub']:
                flag_lines.append(
                    f"# {label} — SPW {spw}, scan {scan}, hub antenna {hub_ant}"
                )
                flag_lines.append(
                    f"mode='manual' scan='{scan}' antenna='{hub_ant}' spw='{spw}'"
                )

    if flag_lines:
        txt_path = f"{vis_base}_spw{spw}_flagcmds.txt"
        with open(txt_path, 'w') as fh:
            fh.write('\n'.join(flag_lines) + '\n')
        print(f"Saved flag commands: {txt_path}")
        print(f"  Apply with: flagdata(vis='{vis}', mode='list', inpfile='{txt_path}')")

    # --- Plots for finite outlier hub antennas ---
    if n_outliers > 0 and output['finite_outliers']:
        # Invert dict: hub_ant -> list of outlier scans
        hub_to_scans = {}
        for scan, info in output['finite_outliers'].items():
            for hub_ant in info['hub']:
                hub_to_scans.setdefault(hub_ant, []).append(scan)

        print(f"\nGenerating plots for {len(hub_to_scans)} hub antenna(s)...")
        for hub_ant, outlier_scans in sorted(hub_to_scans.items()):
            groups = _group_sequential_scans(outlier_scans, all_scans_sorted)
            for grp in groups:
                _plot_hub_weights(vis_base, spw, hub_ant, grp,
                                  all_scans_sorted, ant1_names, ant2_names,
                                  scans, w)
