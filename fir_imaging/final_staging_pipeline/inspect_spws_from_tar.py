"""List the spectral windows of staged MS tars WITHOUT extracting the data.

Only the small SPECTRAL_WINDOW and DATA_DESCRIPTION subtables are pulled out of each tar (the
main table is skipped), so this takes minutes even for a ~300 GB tar and needs almost no scratch.

Use it on the ORIGINAL tars in staged_measurement_sets (job 1 runs combine_spws on the extracted
MS in place, and job 2 would replace the tar with that modified MS), e.g. to see which SPWs
combine_spws dropped as "far-offset" and whether they cover the line.

For every SPW: name, nchan, first/last channel and centre frequency (MHz), channel width (kHz),
number of DATA_DESC rows using it. If a line frequency is known (see below), also the fraction of
the line band (line frequency +/- INSPECT_HALFWIDTH_KMS) that the SPW covers and the channel
index of the line. SPWs are sorted by first-channel frequency; the log of combine_spws numbers
SPWs the same way (SPW ids as in the MS).

Run with (casa does not pass script args reliably, so use env vars):
    INSPECT_TARS="m33_A+B+C+D_oh1665,m33_A+B+C+D_oh1667" \
        casa --nologger --nogui --log2term -c inspect_spws_from_tar.py

Env vars:
    INSPECT_TARS            comma list of tar names (with or without .ms.tar) or full paths.
                            Default: the M33 OH tars (A+B+C+D and C+D, contsub only).
    INSPECT_STAGED_DIR      default /home/ekoch/projects/rrg-eros-ab/ekoch/VLAXL/staged_measurement_sets
    INSPECT_OUTDIR          where <tar>_spws.csv is written (default: current dir)
    INSPECT_LINE_MHZ        override the line frequency (MHz, in the LSRK frame of the MS)
    INSPECT_HALFWIDTH_KMS   half width of the line band, default 100 (M33 rotation)

Default line frequencies (M33, systemic -179 km/s) are the centres of the per-line reference SPW
grids in the M33 listobs: oh1612 1613.199, oh1665 1666.403, oh1667 1668.363, oh1720 1721.564 MHz.
"""

import csv
import os
import re
import shutil
import sys
import tarfile
import tempfile

import numpy as np
from casatools import table

C_KMS = 299792.458
M33_LINE_MHZ = {"oh1612": 1613.1992, "oh1665": 1666.4034, "oh1667": 1668.3634, "oh1720": 1721.5637}
SUBTABLES = ("SPECTRAL_WINDOW", "DATA_DESCRIPTION")

staged_dir = os.environ.get("INSPECT_STAGED_DIR",
                            "/home/ekoch/projects/rrg-eros-ab/ekoch/VLAXL/staged_measurement_sets")
outdir = os.environ.get("INSPECT_OUTDIR", ".")
halfwidth_kms = float(os.environ.get("INSPECT_HALFWIDTH_KMS", "100"))
line_override = os.environ.get("INSPECT_LINE_MHZ", "").strip()
os.makedirs(outdir, exist_ok=True)

tars = [t.strip() for t in os.environ.get("INSPECT_TARS", "").split(",") if t.strip()]
if not tars:
    tars = [f"m33_{c}_{l}" for c in ("A+B+C+D", "C+D") for l in M33_LINE_MHZ]


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


def extract_subtables(path, dest):
    """Extract only SPECTRAL_WINDOW and DATA_DESCRIPTION; return the MS dir name inside dest."""
    ms_dir = None
    with tarfile.open(path, "r") as tf:
        for member in tf:
            parts = member.name.replace("\\", "/").split("/")
            ms_idx = next((i for i, p in enumerate(parts) if p.endswith(".ms")), None)
            if ms_idx is None or len(parts) <= ms_idx + 1 or parts[ms_idx + 1] not in SUBTABLES:
                continue
            ms_dir = parts[ms_idx]
            member.name = "/".join(parts[ms_idx:])
            tf.extract(member, path=dest)
    if ms_dir is None:
        raise RuntimeError(f"{path}: no SPECTRAL_WINDOW subtable found")
    return ms_dir


def band_coverage(f_lo, f_hi, line, halfwidth):
    """Fraction of [line-halfwidth, line+halfwidth] (MHz) inside [f_lo, f_hi]."""
    lo, hi = line - halfwidth, line + halfwidth
    return max(0.0, min(hi, f_hi) - max(lo, f_lo)) / (hi - lo)


def read_spws(ms_path):
    tb = table()
    tb.open(os.path.join(ms_path, "DATA_DESCRIPTION"))
    dd_spw = [int(s) for s in tb.getcol("SPECTRAL_WINDOW_ID")]
    tb.close()
    tb.open(os.path.join(ms_path, "SPECTRAL_WINDOW"))
    out = []
    for i in range(tb.nrows()):
        freqs = np.asarray(tb.getcell("CHAN_FREQ", i)) / 1e6
        widths = np.asarray(tb.getcell("CHAN_WIDTH", i)) / 1e3
        out.append(dict(spw=i, name=str(tb.getcell("NAME", i)), nchan=len(freqs),
                        first_mhz=float(freqs[0]), last_mhz=float(freqs[-1]),
                        ctr_mhz=float(freqs.mean()), width_khz=float(widths[0]),
                        n_dd=dd_spw.count(i), freqs=freqs))
    tb.close()
    return out


for name in tars:
    path = tar_path(name)
    base = os.path.basename(path)[:-len(".ms.tar")]
    print(f"\n=== {base} ===", flush=True)
    work = tempfile.mkdtemp(prefix="inspect_spws_", dir=os.environ.get("SLURM_TMPDIR") or None)
    try:
        ms = os.path.join(work, extract_subtables(path, work))
        spws = read_spws(ms)
    finally:
        shutil.rmtree(work, ignore_errors=True)

    line = line_mhz(path)
    halfwidth = (line * halfwidth_kms / C_KMS) if line else None
    print(f"{len(spws)} SPWs in the SPECTRAL_WINDOW table; "
          + (f"line {line:.4f} MHz +/- {halfwidth:.3f} MHz ({halfwidth_kms:g} km/s)" if line
             else "no line frequency known (coverage columns blank)"))
    print(f"{'spw':>4} {'name':<20} {'nchan':>5} {'first':>10} {'last':>10} {'centre':>10} "
          f"{'wid_kHz':>8} {'nDD':>4} {'line_cov':>8} {'line_ch':>8}")

    rows = []
    for s in sorted(spws, key=lambda s: s["first_mhz"]):
        cov = lchan = None
        if line:
            cov = band_coverage(min(s["first_mhz"], s["last_mhz"]), max(s["first_mhz"], s["last_mhz"]),
                                line, halfwidth)
            if cov > 0:
                lchan = float(np.interp(line, s["freqs"], np.arange(s["nchan"])))
        rows.append(dict(tar=base, spw=s["spw"], name=s["name"], nchan=s["nchan"],
                         first_MHz=round(s["first_mhz"], 4), last_MHz=round(s["last_mhz"], 4),
                         centre_MHz=round(s["ctr_mhz"], 4), width_kHz=round(s["width_khz"], 3),
                         n_datadesc=s["n_dd"],
                         line_band_covered=None if cov is None else round(cov, 3),
                         line_channel=None if lchan is None else round(lchan, 1)))
        print(f"{s['spw']:>4} {s['name']:<20} {s['nchan']:>5} {s['first_mhz']:>10.4f} "
              f"{s['last_mhz']:>10.4f} {s['ctr_mhz']:>10.4f} {s['width_khz']:>8.3f} {s['n_dd']:>4} "
              f"{'' if cov is None else format(cov, '8.2f'):>8} "
              f"{'' if lchan is None else format(lchan, '8.1f'):>8}")

    groups = {}
    for r in rows:
        groups.setdefault((r["name"], r["nchan"], r["width_kHz"]), []).append(r)
    print("By (name, nchan, width): "
          + "; ".join(f"{k[0]} nchan={k[1]} width={k[2]}kHz x{len(v)} "
                      f"(first {min(x['first_MHz'] for x in v):.3f}-{max(x['first_MHz'] for x in v):.3f} MHz)"
                      for k, v in groups.items()))

    csv_path = os.path.join(outdir, f"{base}_spws.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {csv_path}", flush=True)
