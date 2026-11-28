"""Row / flag statistics per (observation, SPW) of a staged MS tar: is the data in each SPW flagged?

Companion to inspect_spws_from_tar.py (which lists the SPW table only). This one reads the main
table, so it extracts the whole tar into a per-job scratch dir (one tar at a time, removed afterwards).
Use the ORIGINAL tar in staged_measurement_sets: job 1's combine_spws has already dropped SPWs from the
MS in its work dir, and job 2 would replace the tar with that modified MS.

Question it answers (M33 17B-162, OH1665/OH1667): combine_spws dropped ~15 "far-offset" SPWs that
cover the line band. Do those SPWs hold unflagged data, and are the rows kept in the reference SPW
for the same observations flagged?

For every (OBSERVATION_ID, SPW) it reports: rows; the fraction of rows with FLAG_ROW set; the fraction
of rows in which every channel/correlation is flagged; the fraction of flagged data points over all
channels and in the line band only (line frequency +/- INSPECT_HALFWIDTH_KMS); the fraction of rows with
infinite WEIGHT. The summary printed at the end groups observations by project and SPWs into the
"reference" SPW (the one with the most rows, as in combine_spws) and "other" SPWs.

Run with (casa does not pass script args reliably, so use env vars; see inspect_spw_rows.sh):
    INSPECT_TARS="m33_A+B+C+D_oh1665" INSPECT_WORKDIR=/path/to/scratch \
        casa --nologger --nogui --log2term -c inspect_spw_rows.py

Env vars:
    INSPECT_TARS            comma list of tar names (with or without .ms.tar) or full paths
                            (default: m33_A+B+C+D_oh1665)
    INSPECT_STAGED_DIR      default /home/ekoch/projects/rrg-eros-ab/ekoch/VLAXL/staged_measurement_sets
    INSPECT_WORKDIR         scratch dir for the extracted MS (default: current dir)
    INSPECT_OUTDIR          where <tar>_spw_rows.csv is written (default: current dir)
    INSPECT_LINE_MHZ        line frequency for the line-band columns (default: M33 OH line from the tar name)
    INSPECT_HALFWIDTH_KMS   half width of the line band, default 100
    INSPECT_KEEP_MS=1       keep the extracted MS (default: remove)
"""

import csv
import os
import re
import shutil
import sys

import numpy as np
from casatools import table

pipe_dir = os.environ.get("INSPECT_PIPE_DIR",
                          "/home/ekoch/lglbs_hi_scripts/fir_imaging/final_staging_pipeline")
sys.path.insert(0, pipe_dir)
from pipeline_common import extract_ms_tar  # noqa: E402

C_KMS = 299792.458
M33_LINE_MHZ = {"oh1612": 1613.1992, "oh1665": 1666.4034, "oh1667": 1668.3634, "oh1720": 1721.5637}
PROJECT_LABELS = {"uid://evla/pdb/29631320": "14B-088", "uid://evla/pdb/34051870": "17B-162",
                  "uid://evla/pdb/38078528": "20A-346"}     # M33 (as logged in listobs)
flag_chunk_elems = 2e8     # max FLAG elements read at once

staged_dir = os.environ.get("INSPECT_STAGED_DIR",
                            "/home/ekoch/projects/rrg-eros-ab/ekoch/VLAXL/staged_measurement_sets")
workdir = os.environ.get("INSPECT_WORKDIR", ".")
outdir = os.environ.get("INSPECT_OUTDIR", ".")
halfwidth_kms = float(os.environ.get("INSPECT_HALFWIDTH_KMS", "100"))
line_override = os.environ.get("INSPECT_LINE_MHZ", "").strip()
keep_ms = os.environ.get("INSPECT_KEEP_MS", "0").strip().lower() in ("1", "true", "yes")
os.makedirs(workdir, exist_ok=True)
os.makedirs(outdir, exist_ok=True)
tars = [t.strip() for t in os.environ.get("INSPECT_TARS", "m33_A+B+C+D_oh1665").split(",") if t.strip()]


def tar_path(name):
    if os.path.isabs(name):
        return name
    for cand in (name, name + ".tar", name + ".ms.tar"):
        p = os.path.join(staged_dir, cand)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"{name}: no tar found in {staged_dir}")


def line_mhz(path):
    if line_override:
        return float(line_override)
    m = re.match(r"m33_.*_(oh\d+)(_nocontsub)?\.ms\.tar$", os.path.basename(path))
    return M33_LINE_MHZ.get(m.group(1)) if m else None


def accumulate(acc, obs, frow, flag, weight, ch_slice):
    """Add one chunk to acc = {obs: np.array([rows, frow, allflag, flagged_pts, pts, lb_flagged, lb_pts, inf])}.

    obs (n,), frow (n,) bool, flag (ncorr, nchan, n) bool, weight (ncorr, n) or None,
    ch_slice = slice of the line-band channels (or None).
    """
    n = len(obs)
    ncorr, nchan, _ = flag.shape
    row_flagged = flag.sum(axis=(0, 1))                                   # flagged points per row
    all_flag = row_flagged == ncorr * nchan
    if ch_slice is not None:
        lb_flagged = flag[:, ch_slice, :].sum(axis=(0, 1))
        lb_pts = ncorr * (ch_slice.stop - ch_slice.start)
    else:
        lb_flagged, lb_pts = np.zeros(n, dtype=np.int64), 0
    inf_w = np.isinf(weight).any(axis=0) if weight is not None else np.zeros(n, dtype=bool)
    uniq, inv = np.unique(obs, return_inverse=True)
    stats = np.vstack([np.ones(n), frow, all_flag, row_flagged, np.full(n, ncorr * nchan),
                       lb_flagged, np.full(n, lb_pts), inf_w])
    for k in range(len(uniq)):
        acc.setdefault(int(uniq[k]), np.zeros(8))
        acc[int(uniq[k])] += stats[:, inv == k].sum(axis=1)


def scan_dd(tb, ms, dd, nchan, ch_slice):
    acc = {}
    sub = tb.query(f"DATA_DESC_ID=={dd}")
    try:
        nrow = sub.nrows()
        if nrow == 0:
            return acc
        ncorr = np.shape(sub.getcell("FLAG", 0))[0]
        step = max(1000, int(flag_chunk_elems // (ncorr * nchan)))
        has_w = "WEIGHT" in sub.colnames()
        for start in range(0, nrow, step):
            n = min(step, nrow - start)
            obs = sub.getcol("OBSERVATION_ID", startrow=start, nrow=n)
            frow = sub.getcol("FLAG_ROW", startrow=start, nrow=n)
            flag = sub.getcol("FLAG", startrow=start, nrow=n)              # (ncorr, nchan, n)
            weight = sub.getcol("WEIGHT", startrow=start, nrow=n) if has_w else None
            accumulate(acc, obs, frow, flag, weight, ch_slice)
    finally:
        sub.close()
    return acc


def pct(a, b):
    return 100.0 * a / b if b else float("nan")


for name in tars:
    path = tar_path(name)
    base = os.path.basename(path)[:-len(".ms.tar")]
    print(f"\n=== {base} ===", flush=True)
    line = line_mhz(path)
    halfwidth = line * halfwidth_kms / C_KMS if line else None
    print(f"Extracting {path} -> {workdir}", flush=True)
    ms = extract_ms_tar(path, workdir)
    try:
        tb = table()
        tb.open(os.path.join(ms, "OBSERVATION"))
        project = [str(p) for p in tb.getcol("PROJECT")]
        tb.close()
        tb.open(os.path.join(ms, "DATA_DESCRIPTION"))
        dd_spw = [int(s) for s in tb.getcol("SPECTRAL_WINDOW_ID")]
        tb.close()
        tb.open(os.path.join(ms, "SPECTRAL_WINDOW"))
        spw_info = {}
        for i in range(tb.nrows()):
            f = np.asarray(tb.getcell("CHAN_FREQ", i)) / 1e6
            spw_info[i] = dict(nchan=len(f), first=float(f[0]), last=float(f[-1]), freqs=f)
        tb.close()

        tb.open(ms)
        results = {}                                   # (obs, spw) -> acc row
        for dd, spw in enumerate(dd_spw):
            info = spw_info[spw]
            ch_slice = None
            if line:
                idx = np.where(np.abs(info["freqs"] - line) <= halfwidth)[0]
                if idx.size:
                    ch_slice = slice(int(idx.min()), int(idx.max()) + 1)
            print(f"  DD {dd} (SPW {spw}, {info['first']:.3f}-{info['last']:.3f} MHz): "
                  f"line band channels {ch_slice.start}-{ch_slice.stop - 1}" if ch_slice else
                  f"  DD {dd} (SPW {spw}, {info['first']:.3f}-{info['last']:.3f} MHz): "
                  f"{'line band not covered' if line else 'no line frequency'}", flush=True)
            for o, a in scan_dd(tb, ms, dd, info["nchan"], ch_slice).items():
                results.setdefault((o, spw), np.zeros(8))
                results[(o, spw)] += a
        tb.close()
    finally:
        if not keep_ms and os.path.exists(ms):
            shutil.rmtree(ms)

    rows_per_spw = {}
    for (o, spw), a in results.items():
        rows_per_spw[spw] = rows_per_spw.get(spw, 0) + a[0]
    ref = max(rows_per_spw, key=rows_per_spw.get)
    label = lambda o: PROJECT_LABELS.get(project[o] if o < len(project) else "?", project[o] if o < len(project) else "?")

    csv_path = os.path.join(outdir, f"{base}_spw_rows.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["obs", "project", "spw", "spw_role", "spw_first_MHz", "spw_last_MHz", "rows",
                    "pct_flag_row", "pct_rows_all_flagged", "pct_points_flagged",
                    "pct_lineband_flagged", "pct_inf_weight"])
        for (o, spw) in sorted(results):
            a = results[(o, spw)]
            w.writerow([o, label(o), spw, "ref" if spw == ref else "other",
                        round(spw_info[spw]["first"], 4), round(spw_info[spw]["last"], 4), int(a[0]),
                        round(pct(a[1], a[0]), 2), round(pct(a[2], a[0]), 2), round(pct(a[3], a[4]), 2),
                        round(pct(a[5], a[6]), 2) if a[6] else "", round(pct(a[7], a[0]), 2)])
    print(f"Wrote {csv_path}")

    # ---- summary by project x (ref | other)
    grp = {}
    for (o, spw), a in results.items():
        g = grp.setdefault((label(o), "ref" if spw == ref else "other"), np.zeros(8))
        g += a
    print(f"\nReference SPW (most rows) = SPW {ref}. Totals by project and SPW role:")
    print(f"{'project':<10} {'role':<6} {'rows':>11} {'%rowflag':>9} {'%allflag':>9} {'%pts':>7} "
          f"{'%lineband':>9} {'%infW':>7}")
    for (p, role), a in sorted(grp.items()):
        print(f"{p:<10} {role:<6} {int(a[0]):>11d} {pct(a[1], a[0]):>9.1f} {pct(a[2], a[0]):>9.1f} "
              f"{pct(a[3], a[4]):>7.1f} {(format(pct(a[5], a[6]), '9.1f') if a[6] else '        -'):>9} "
              f"{pct(a[7], a[0]):>7.1f}")
    print("(%allflag = rows with every channel flagged; %lineband = flagged points in the line band only)")
