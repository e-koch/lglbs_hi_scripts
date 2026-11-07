"""Flagged percentage vs channel for a CASA MS, one plot (+ CSV) per SPW in use.

The percentage is the fraction of all data points (rows x correlations) in a channel that
are flagged, so it matches the totals in <base>_flag_summary.txt. Run it after the flags
have been applied (find_weight_outliers.py) to see the final state of the staged MS.

Each PNG has two panels sharing the channel axis (channel number on the bottom axis, channel
frequency in MHz on the top axis): linear 0-100% on top, log scale below so that small
fractions and isolated flagged channels are both visible.

Run with:
    casa --nologger --nogui --log2term -c plot_flag_vs_channel.py <vis> [outdir]

Outputs go in [outdir] if given (created if needed), else next to the MS:
    <base>_spw<N>_flag_vs_channel.png
    <base>_spw<N>_flag_vs_channel.csv    columns: spw, chan, freq_MHz, flagged_pct

Frequencies are the channel centres in the MS frame (LSRK for the staged MSs).
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np
import os
import sys

from casatools import msmetadata, table

flag_chunk_elems = 2e8     # max FLAG elements (rows x corr x chan) read at once
attention_pct = 50.0       # channels flagged at least this much are listed in the log


def script_args(script_name):
    """Arguments following `script_name` on a `casa -c script args` command line."""
    idx = max(i for i, a in enumerate(sys.argv) if a.endswith(script_name))
    return sys.argv[idx + 1:]


def flagged_per_channel(tb, nchan):
    """(flagged data points per channel, total data points per channel) for an open table."""
    nrow = tb.nrows()
    ncorr = np.shape(tb.getcell("FLAG", 0))[0]
    step = max(1000, int(flag_chunk_elems // (ncorr * nchan)))
    counts = np.zeros(nchan, dtype=np.int64)
    for start in range(0, nrow, step):
        n = min(step, nrow - start)
        flag = tb.getcol("FLAG", startrow=start, nrow=n)      # (ncorr, nchan, n)
        counts += flag.sum(axis=(0, 2))
    return counts, nrow * ncorr


def make_plot(png, title, chan, freq_mhz, pct, floored, ylow, frame):
    fig, (ax_lin, ax_log) = plt.subplots(2, 1, figsize=(11, 6.5), sharex=True,
                                         gridspec_kw=dict(height_ratios=[1, 1], hspace=0.08))
    for ax in (ax_lin, ax_log):
        ax.step(chan, pct if ax is ax_lin else floored, where="mid", color="C0", lw=1.0)
        ax.grid(alpha=0.3)
    ax_lin.set_ylim(-2, 102)
    ax_lin.set_ylabel("flagged (%)")
    ax_log.set_yscale("log")
    ax_log.set_ylim(ylow, 130)
    ax_log.set_ylabel("flagged (%), log")
    ax_log.set_xlabel("channel number")
    ax_log.set_xlim(chan[0] - 0.5, chan[-1] + 0.5)
    if (pct == 0).any():    # zero-flagged channels are drawn at the bottom of the log panel
        ax_log.text(0.995, 0.03, f"zero plotted at {ylow:.2g}%", transform=ax_log.transAxes,
                    ha="right", va="bottom", fontsize=8, color="0.4")

    if len(chan) > 1:
        df = (freq_mhz[-1] - freq_mhz[0]) / (len(chan) - 1)    # MHz per channel
        f0 = freq_mhz[0]
        top = ax_lin.secondary_xaxis("top", functions=(lambda c: f0 + c * df,
                                                       lambda f: (f - f0) / df))
        top.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))   # full MHz, no "+1.42e3"
        top.set_xlabel(f"frequency (MHz{', ' + frame if frame else ''})")
    ax_lin.set_title(title, fontsize=10)
    fig.savefig(png, dpi=110, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    args = script_args("plot_flag_vs_channel.py")
    vis = args[0].rstrip("/").rstrip(os.sep)
    outdir = args[1] if len(args) > 1 else os.path.dirname(os.path.abspath(vis))
    os.makedirs(outdir, exist_ok=True)
    base = os.path.join(outdir, os.path.splitext(os.path.basename(vis))[0])

    tb = table()
    tb.open(vis + "/DATA_DESCRIPTION")
    dd_spw = tb.getcol("SPECTRAL_WINDOW_ID")
    tb.close()
    tb.open(vis)
    dd_rows = tb.getcol("DATA_DESC_ID")
    tb.close()
    dds_in_use = [int(d) for d in np.unique(dd_rows)]
    spws = sorted({int(dd_spw[d]) for d in dds_in_use})

    msmd = msmetadata()
    msmd.open(vis)
    for spw in spws:
        nchan = int(msmd.nchan(spw))
        freq_mhz = np.asarray(msmd.chanfreqs(spw)) / 1e6
        try:
            frame = str(msmd.measfreqframe(spw))
        except Exception:
            frame = ""

        flagged = np.zeros(nchan, dtype=np.int64)
        total = 0
        tb.open(vis)
        for dd in (d for d in dds_in_use if int(dd_spw[d]) == spw):
            sub = tb.query(f"DATA_DESC_ID=={dd}")
            c, t = flagged_per_channel(sub, nchan)
            sub.close()
            flagged += c
            total += t
        tb.close()

        pct = 100.0 * flagged / total
        chan = np.arange(nchan)

        # log-panel data: zeros are drawn at a floor just below the smallest non-zero value
        positive = pct[pct > 0]
        ylow = max(1e-4, 0.3 * positive.min()) if positive.size else 1e-3
        floored = np.maximum(pct, ylow)

        tag = f"{base}_spw{spw}_flag_vs_channel"
        title = (f"{os.path.basename(vis)}  SPW {spw}: {nchan} channels, "
                 f"overall {100.0 * flagged.sum() / (total * nchan):.3f}% flagged, "
                 f"median channel {np.median(pct):.3f}%, max {pct.max():.2f}% (ch {int(pct.argmax())})")
        make_plot(tag + ".png", title, chan, freq_mhz, pct, floored, ylow, frame)
        np.savetxt(tag + ".csv", np.c_[np.full(nchan, spw), chan, freq_mhz, pct],
                   delimiter=",", header="spw,chan,freq_MHz,flagged_pct", comments="",
                   fmt=["%d", "%d", "%.6f", "%.4f"])
        print(f"SPW {spw}: {title.split(': ', 1)[1]}")
        bad = np.where(pct >= attention_pct)[0]
        if bad.size:
            print(f"  channels flagged >= {attention_pct:g}%: "
                  + ", ".join(f"ch{c} ({pct[c]:.1f}%, {freq_mhz[c]:.3f} MHz)" for c in bad[:20])
                  + (" ..." if bad.size > 20 else ""))
        print(f"Wrote {tag}.png and .csv")
    msmd.close()
