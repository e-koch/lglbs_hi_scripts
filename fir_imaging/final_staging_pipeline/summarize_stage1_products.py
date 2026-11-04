#!/usr/bin/env python
"""Summary reports + PNG galleries for the stage-1 (job 1) products of the final staging pipeline.

Job 1 (job1_prep_and_test.sh) leaves one ``<name>_products/`` directory per staged MS:

    <name>_flag_summary.txt            from find_weight_outliers.py (flagging before/after)
    <name>[_spwN]_flagcmds.txt         hub-antenna flag commands that were applied
    <name>.listobs.txt                 listobs of the MS
    test_images/test_imaging_summary.txt
    test_images/test_<field>_<cfg>.{psf,residual}      CASA images from run_test_imaging.py

For every products directory this script writes (to ``--outdir``):

    <name>_summary.md      flagging + imaging summary and automatic attention flags
    <name>_gallery.png     PSF + residual for every field x VLA config
and, across all runs:

    overview.md            one table row per run plus the list of runs needing attention
    overview.csv           the same numbers in machine-readable form

Only the standard scientific Python stack is needed (numpy, matplotlib). CASA images are
read with casatools if it is importable (i.e. inside CASA), otherwise with casa-formats-io.

Usage:
    python summarize_stage1_products.py [DATA_DIR] [-o OUTDIR] [--match GLOB] [--no-gallery]
        DATA_DIR  directory holding the *_products dirs (default /home/koche/data)
        -o        output directory (default DATA_DIR/stage1_qa_reports)
        --match   only process products dirs whose name matches this glob, e.g. 'ic1613*oh*'
        --no-gallery  skip the PNGs (text reports only)
"""

import argparse
import csv
import fnmatch
import glob
import os
import re
import sys
import types
from collections import OrderedDict, defaultdict

import numpy as np

DEFAULT_DATA_DIR = "/home/koche/data"

# --- thresholds for the automatic "attention" flags -------------------------------------
HIGH_FLAG_PCT = 2.0         # total flagged fraction after hub-antenna flagging (%)
BIG_FLAG_CHANGE_PCT = 0.2   # flagging added by the hub-antenna step (%)
FULLY_FLAGGED_SCAN_PCT = 99.0   # a scan this flagged after flagging is effectively lost (%)
ANT_FLAG_FACTOR = 5.0       # antenna flagged > this x the median antenna ...
ANT_FLAG_MIN_PCT = 1.0      # ... and above this absolute level (%)
RMS_OUTLIER_FACTOR = 1.3    # residual rms vs median of the same config
PSF_PEAK_TOL = 1e-3         # |psf_peak - 1| tolerance
RESID_PEAK_SNR_MAX = 8.0    # |residual| peak / rms; much more than this is not noise-like
PSF_DISPLAY = 0.2           # gallery PSF stretch: +-this
RESID_DISPLAY_SIGMA = 4.0   # gallery residual stretch: +-this x rms


# =========================================================================================
# CASA image reading
# =========================================================================================

_image_reader = None


def _make_reader():
    """Return f(path) -> 2D float array for a single-plane CASA image."""
    try:  # inside CASA, or any env with casatools
        from casatools import image as iatool

        def read_casatools(path):
            ia = iatool()
            ia.open(path)
            try:
                return np.squeeze(ia.getchunk()).T.astype(float)  # (x, y) -> (y, x)
            finally:
                ia.close()
        return read_casatools
    except ImportError:
        pass

    try:
        from casa_formats_io import image_to_dask
    except ImportError as exc:
        # casa-formats-io has a compiled helper (_casa_chunking); a binary built for another
        # OS/Python raises ImportError here. Substitute a numpy version of it (checked
        # bit-for-bit against the compiled one for float32/float64 data).
        if "_casa_chunking" in str(exc) or "ELF" in str(exc) or "invalid" in str(exc):
            stub = types.ModuleType("casa_formats_io._casa_chunking")

            def _combine_chunks(array_1d, itemsize, n4, n3, n2, n1, o4, o3, o2, o1):
                """Interleave `oversample` native tiles (stored one after another, each
                C-ordered over (x, y, stokes, chan)) into one array C-ordered the same way."""
                nat, ov = (n1, n2, n3, n4), (o1, o2, o3, o4)
                if all(o == 1 for o in ov):
                    return array_1d
                if itemsize == 1:
                    raise NotImplementedError("bit-packed masks need the compiled helper")
                tiles = np.frombuffer(array_1d, dtype=np.uint8).reshape(ov + nat + (itemsize,))
                full = tiles.transpose(0, 4, 1, 5, 2, 6, 3, 7, 8)
                return np.ascontiguousarray(full).reshape(-1)
            stub._combine_chunks = _combine_chunks
            sys.modules["casa_formats_io._casa_chunking"] = stub
            from casa_formats_io import image_to_dask
        else:
            raise ImportError("Need casatools or casa-formats-io to read CASA images") from exc

    def read_cfio(path):
        return np.asarray(image_to_dask(path)).squeeze().astype(float)
    return read_cfio


def read_image(path):
    global _image_reader
    if _image_reader is None:
        _image_reader = _make_reader()
    return _image_reader(path)


# =========================================================================================
# Parsers
# =========================================================================================

def parse_flag_summary(path):
    """Parse *_flag_summary.txt -> dict (None if the file is missing)."""
    if not os.path.exists(path):
        return None
    out = dict(ms=None, observation=None, cmd_files=[], before=None, after=None, change=None,
               npts=None, sections={})
    section = None
    with open(path) as fh:
        lines = fh.read().splitlines()
    in_cmds = False
    for line in lines:
        s = line.strip()
        if line.startswith("Flag summary for"):
            out["ms"] = line.split("for", 1)[1].strip()
        elif line.startswith("Observation:"):
            out["observation"] = line.split(":", 1)[1].strip()
        elif line.startswith("Flag command files applied"):
            in_cmds = True
        elif line.startswith("Total flagged (before)"):
            in_cmds = False
            m = re.search(r":\s*([\d.]+)\s*%\s*\((\d+)\s*/\s*(\d+)", line)
            out["before"], out["npts"] = float(m.group(1)), int(m.group(3))
        elif line.startswith("Total flagged (after)"):
            out["after"] = float(re.search(r":\s*([\d.]+)\s*%", line).group(1))
        elif line.startswith("Change"):
            out["change"] = float(re.search(r":\s*([+-]?[\d.]+)\s*%", line).group(1))
        elif line.startswith("Per "):
            section = line.split("Per ", 1)[1].strip()
            out["sections"][section] = OrderedDict()
        elif in_cmds and s and not s.startswith("(none"):
            out["cmd_files"].append(s)
        elif section and s and not s.startswith("name"):
            parts = s.split()
            if len(parts) == 4:
                try:
                    out["sections"][section][parts[0]] = tuple(float(x) for x in parts[1:])
                except ValueError:
                    pass
    return out


def parse_flagcmds(paths):
    """Summarise hub-antenna flag command files -> dict of counts/sets."""
    cmds, scans, ants, spws = 0, set(), set(), set()
    for p in paths:
        with open(p) as fh:
            for line in fh:
                if line.startswith("mode="):
                    cmds += 1
                    for key, store in (("scan", scans), ("antenna", ants), ("spw", spws)):
                        m = re.search(key + r"='([^']*)'", line)
                        if m:
                            store.update(x for x in m.group(1).split(",") if x)
    return dict(n_cmds=cmds, scans=scans, antennas=ants, spws=spws)


def parse_imaging_summary(path):
    """Parse test_imaging_summary.txt -> dict (None if the file is missing)."""
    if not os.path.exists(path):
        return None
    out = dict(ms=None, spw=None, spws=[], chan=None, nchan=None, freq_mhz=None, datacolumn=None,
               configs=OrderedDict(), rows=[])
    with open(path) as fh:
        lines = fh.read().splitlines()
    for line in lines:
        m = re.match(r"Test imaging for (.*)", line)
        if m:
            out["ms"] = m.group(1).strip()
        # "SPW 0 channel ..." (first-SPW-only imaging) or "SPW 0,1 channel ..." (all SPWs)
        m = re.match(r"SPWs? ([\d,]+) channel (\d+)/(\d+) \(([\d.]+) MHz\).*datacolumn=(\w+)", line)
        if m:
            out["spw"] = m.group(1)                       # label, e.g. "0" or "0,1"
            out["spws"] = [int(x) for x in m.group(1).split(",")]
            out["chan"], out["nchan"] = int(m.group(2)), int(m.group(3))
            out["freq_mhz"], out["datacolumn"] = float(m.group(4)), m.group(5)
        m = re.match(r"^\s*([ABCD])\s+(\S+)\s+(\d+)\s*$", line)
        if m:
            out["configs"][m.group(1)] = dict(obsids=[int(x) for x in m.group(2).split(",")],
                                              p95=int(m.group(3)))
        m = re.match(r"^(\S+)\s+([ABCD])\s+([\d.]+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(ok|FAILED.*)$", line)
        if m:
            out["rows"].append(dict(field=m.group(1), config=m.group(2), cell=float(m.group(3)),
                                    imsize=int(m.group(4)), psf_peak=float(m.group(5)),
                                    resid_rms=float(m.group(6)), status=m.group(7)))
    return out


def parse_listobs(path):
    """Parse listobs -> {(obsid, field_name): set(spw)} plus field ids / names.

    Returns dict(rows={(obsid, field_name): nrows_by_spw dict}, fields={name: id}).
    Field names come from the Fields table (listobs truncates them in the scan table, so
    scans are matched on FldId).
    """
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        lines = fh.read().splitlines()
    fields = OrderedDict()
    in_fields = False
    for line in lines:
        if line.startswith("Fields:"):
            in_fields = True
            continue
        if in_fields:
            if line.startswith("Spectral Windows"):
                break
            m = re.match(r"^\s+(\d+)\s+\S+\s+(\S+)\s+\d\d:", line)
            if m:
                fields[m.group(2)] = int(m.group(1))
    id_to_name = {v: k for k, v in fields.items()}
    rows = defaultdict(lambda: defaultdict(int))   # (obsid, name) -> spw -> nrows
    obsid = None
    for line in lines:
        m = re.match(r"^\s+ObservationID = (\d+)", line)
        if m:
            obsid = int(m.group(1))
            continue
        m = re.search(r"\d\.\d\s+(\d+)\s+(\d+)\s+(.+?)\s{2,}(\d+)\s+\[([\d,\s]+)\]", line)
        if m and obsid is not None:
            fld = int(m.group(2))
            nrows = int(m.group(4))
            spws = [int(x) for x in m.group(5).split(",")]
            for spw in spws:
                rows[(obsid, id_to_name.get(fld, m.group(3).strip()))][spw] += nrows // len(spws)
    return dict(rows=rows, fields=fields)


# =========================================================================================
# Per-run analysis
# =========================================================================================

def run_label(products_dir):
    return os.path.basename(products_dir.rstrip("/"))[: -len("_products")]


def explain_failure(row, imaging, listobs):
    """Classify a FAILED row -> (kind, text). kind is 'expected' or 'unexplained'."""
    if listobs is None:
        return "unexplained", "no listobs to check against"
    cfg_obsids = imaging["configs"].get(row["config"], {}).get("obsids", [])
    spws = imaging["spws"]
    spw = imaging["spw"]
    have = defaultdict(int)       # spw -> rows for this field in this config's obsids
    for (oid, fname), by_spw in listobs["rows"].items():
        if fname == row["field"] and oid in cfg_obsids:
            for s, n in by_spw.items():
                have[s] += n
    if any(have.get(sp, 0) > 0 for sp in spws):
        return "unexplained", f"listobs shows data in imaged SPW {spw} for this field/config"
    if not have:
        return "expected", f"field not observed in {row['config']}-config obsids"
    others = ", ".join(str(s) for s in sorted(have))
    return "expected", f"field only observed in SPW {others}, test image uses SPW {spw}"


def image_stats(path_resid, path_psf):
    """Pixel statistics read from the CASA images."""
    out = {}
    r = read_image(path_resid)
    good = np.isfinite(r) & (r != 0)
    rms = float(np.sqrt(np.mean(r[good] ** 2))) if good.any() else np.nan
    out["rms"] = rms
    out["resid_peak_snr"] = float(np.max(np.abs(r[good])) / rms) if good.any() else np.nan
    out["resid_mean_over_rms"] = float(np.mean(r[good]) / rms) if good.any() else np.nan
    p = read_image(path_psf)
    out["psf_max"] = float(np.nanmax(p))
    out["psf_min"] = float(np.nanmin(p))
    return out


def analyse_run(products_dir):
    """Collect everything about one products dir."""
    label = run_label(products_dir)
    run = dict(label=label, dir=products_dir)
    parts = label.split("_")
    run["target"] = parts[0]
    run["nocontsub"] = label.endswith("_nocontsub")
    run["line"] = "_".join(p for p in parts[2:] if p != "nocontsub")
    run["flag"] = parse_flag_summary(os.path.join(products_dir, f"{label}_flag_summary.txt"))
    cmd_files = sorted(glob.glob(os.path.join(products_dir, "*_flagcmds.txt")))
    # <name>_flagcmds.txt is the union of the per-SPW <name>_spwN_flagcmds.txt files
    union = [f for f in cmd_files if "_spw" not in os.path.basename(f)]
    run["flagcmd_files"] = [os.path.basename(f) for f in cmd_files]
    run["flagcmds"] = parse_flagcmds(union or cmd_files)
    img_dir = os.path.join(products_dir, "test_images")
    run["img_dir"] = img_dir
    run["imaging"] = parse_imaging_summary(os.path.join(img_dir, "test_imaging_summary.txt"))
    run["listobs"] = parse_listobs(os.path.join(products_dir, f"{label}.listobs.txt"))

    attention = []
    flag = run["flag"]
    if flag is None:
        attention.append("MISSING flag summary")
    else:
        if flag["after"] is not None and flag["after"] > HIGH_FLAG_PCT:
            attention.append(f"High flagged fraction after flagging: {flag['after']:.2f}%")
        if flag["change"] is not None and flag["change"] > BIG_FLAG_CHANGE_PCT:
            attention.append(f"Hub-antenna flagging added {flag['change']:.3f}% "
                             f"({flag['before']:.3f}% -> {flag['after']:.3f}%)")
        scans = flag["sections"].get("scan", {})
        dead = [k for k, v in scans.items() if v[1] >= FULLY_FLAGGED_SCAN_PCT]
        if dead:
            attention.append(f"{len(dead)} scan(s) >= {FULLY_FLAGGED_SCAN_PCT:g}% flagged after flagging: "
                             + ", ".join(f"{k} ({scans[k][0]:.1f}% -> {scans[k][1]:.1f}%)" for k in dead[:10])
                             + (" ..." if len(dead) > 10 else ""))
        ants = flag["sections"].get("antenna", {})
        if len(ants) > 3:
            med = float(np.median([v[1] for v in ants.values()]))
            heavy = [(k, v[1]) for k, v in ants.items()
                     if v[1] > max(ANT_FLAG_MIN_PCT, ANT_FLAG_FACTOR * med)]
            if heavy:
                attention.append(f"Antenna(s) heavily flagged vs median {med:.3f}%: "
                                 + ", ".join(f"{k} ({v:.2f}%)" for k, v in sorted(heavy, key=lambda x: -x[1])))

    imaging = run["imaging"]
    run["n_ok"] = run["n_fail_expected"] = run["n_fail_unexplained"] = 0
    run["median_rms"] = {}
    if imaging is None:
        attention.append("MISSING test_imaging_summary.txt")
    else:
        for row in imaging["rows"]:
            row["stem"] = f"test_{row['field']}_{row['config']}"
            row["psf_path"] = os.path.join(img_dir, row["stem"] + ".psf")
            row["resid_path"] = os.path.join(img_dir, row["stem"] + ".residual")
            row["ok"] = row["status"] == "ok"
            row["images_exist"] = os.path.isdir(row["psf_path"]) and os.path.isdir(row["resid_path"])
            row["note"] = ""
            if row["ok"]:
                run["n_ok"] += 1
                if not row["images_exist"]:
                    row["note"] = "image files missing on disk"
                    attention.append(f"{row['field']} {row['config']}: reported ok but image files missing")
                else:
                    try:
                        row.update(image_stats(row["resid_path"], row["psf_path"]))
                    except Exception as exc:   # one unreadable image must not stop the report
                        row["note"] = f"could not read images: {exc}"
                        attention.append(f"{row['field']} {row['config']}: {row['note']}")
            else:
                kind, why = explain_failure(row, imaging, run["listobs"])
                row["fail_kind"], row["note"] = kind, why
                run["n_fail_" + kind] += 1
                if kind == "unexplained":
                    attention.append(f"{row['field']} {row['config']}: FAILED, {why}")
        # per-config residual rms outliers, PSF / residual sanity (one summary line per kind)
        rms_out, psf_bad, peaky = [], [], []
        for cfg in imaging["configs"]:
            oks = [r for r in imaging["rows"] if r["config"] == cfg and r["ok"]
                   and np.isfinite(r.get("rms", np.nan))]
            if oks:
                med = float(np.median([r["rms"] for r in oks]))
                run["median_rms"][cfg] = med
                rms_out += [f"{r['field']} {cfg} ({r['rms'] * 1e3:.2f} mJy/beam, {r['rms'] / med:.2f}x)"
                            for r in oks if r["rms"] > RMS_OUTLIER_FACTOR * med]
        for r in imaging["rows"]:
            if r["ok"] and "psf_max" in r:
                if abs(r["psf_max"] - 1.0) > PSF_PEAK_TOL:
                    psf_bad.append(f"{r['field']} {r['config']} ({r['psf_max']:.4f})")
                if r["resid_peak_snr"] > RESID_PEAK_SNR_MAX:
                    peaky.append((r["resid_peak_snr"], f"{r['field']} {r['config']}"))
        if rms_out:
            attention.append(f"Residual rms > {RMS_OUTLIER_FACTOR}x the per-config median: " + "; ".join(rms_out))
        if psf_bad:
            attention.append("PSF peak != 1: " + "; ".join(psf_bad))
        if peaky:
            peaky.sort(reverse=True)
            attention.append(
                f"{len(peaky)} image(s) have a residual peak > {RESID_PEAK_SNR_MAX:g} sigma (bright "
                f"emission/continuum, RFI or a bad scan; expected for nocontsub runs): "
                + ", ".join(f"{n} ({v:.0f} sigma)" for v, n in peaky[:6])
                + (" ..." if len(peaky) > 6 else ""))
        n_unexp = run["n_fail_unexplained"]
        if run["n_fail_expected"]:
            attention.append(f"{run['n_fail_expected']} of {len(imaging['rows'])} field/config test images could not "
                             f"be made: no data in the imaged SPW/config (see table); coverage is incomplete")
        if run["n_ok"] == 0 and imaging["rows"]:
            attention.append("NO test images were produced")
    run["attention"] = attention
    return run


# =========================================================================================
# Gallery
# =========================================================================================

def make_gallery(run, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    imaging = run["imaging"]
    if imaging is None or not imaging["rows"]:
        return False
    fields = list(OrderedDict.fromkeys(r["field"] for r in imaging["rows"]))
    configs = sorted(imaging["configs"]) if imaging["configs"] else sorted({r["config"] for r in imaging["rows"]})
    lookup = {(r["field"], r["config"]): r for r in imaging["rows"]}
    ncols = 2 * len(configs)
    cell = 2.6
    fig, axes = plt.subplots(len(fields), ncols, figsize=(cell * ncols + 1.6, cell * len(fields) + 1.2),
                             squeeze=False)
    for i, field in enumerate(fields):
        for j, cfg in enumerate(configs):
            row = lookup.get((field, cfg))
            for k, kind in enumerate(("psf", "residual")):
                ax = axes[i][2 * j + k]
                ax.set_xticks([])
                ax.set_yticks([])
                if i == 0:
                    ax.set_title(f"{cfg}-config  {'PSF' if kind == 'psf' else 'residual'}", fontsize=10)
                if j == 0 and k == 0:
                    ax.set_ylabel(field, fontsize=9)
                if row is None or not row["ok"] or "rms" not in row:
                    msg = "not imaged" if row is None else ("FAILED\n(no data)" if not row["ok"] else "unreadable")
                    ax.text(0.5, 0.5, msg, ha="center", va="center", transform=ax.transAxes,
                            fontsize=9, color="firebrick")
                    for sp in ax.spines.values():
                        sp.set_color("firebrick")
                    continue
                img = read_image(row["psf_path"] if kind == "psf" else row["resid_path"])
                if kind == "psf":
                    ax.imshow(img, origin="lower", cmap="gray", vmin=-PSF_DISPLAY, vmax=PSF_DISPLAY)
                    ax.text(0.03, 0.03, f"{row['imsize']}px, {row['cell']:.1f}\"/px", color="yellow",
                            fontsize=7, transform=ax.transAxes)
                else:
                    lim = RESID_DISPLAY_SIGMA * row["rms"]
                    ax.imshow(img, origin="lower", cmap="RdBu_r", vmin=-lim, vmax=lim)
                    ax.text(0.03, 0.03, f"rms {row['rms'] * 1e3:.2f} mJy/bm", color="black", fontsize=7,
                            transform=ax.transAxes,
                            bbox=dict(facecolor="white", alpha=0.7, pad=1, lw=0))
    flag = run["flag"]
    flagtxt = (f"flagged {flag['before']:.3f}% -> {flag['after']:.3f}%" if flag and flag["after"] is not None
               else "no flag summary")
    spwtxt = (f"SPW {imaging['spw']} ch {imaging['chan']}/{imaging['nchan']} ({imaging['freq_mhz']:.3f} MHz)"
              if imaging["spw"] is not None else "")
    fig.suptitle(f"{run['label']}   {spwtxt}   {flagtxt}\n"
                 f"PSF stretch +-{PSF_DISPLAY}; residual stretch +-{RESID_DISPLAY_SIGMA:g} sigma (per image)",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.985 if len(fields) > 4 else 0.96))
    fig.savefig(out_png, dpi=80)
    plt.close(fig)
    return True


# =========================================================================================
# Report writers
# =========================================================================================

def md_table(headers, rows, align=None):
    align = align or ["l"] * len(headers)
    sep = [":--" if a == "l" else "--:" for a in align]
    out = ["| " + " | ".join(headers) + " |", "| " + " | ".join(sep) + " |"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def fmt_list(items, maxn=12):
    items = sorted(items, key=lambda x: (not str(x).isdigit(), int(x) if str(x).isdigit() else str(x)))
    s = ", ".join(str(i) for i in items[:maxn])
    return s + (f", ... ({len(items)} total)" if len(items) > maxn else "")


def write_run_report(run, path, gallery_name):
    flag, imaging = run["flag"], run["imaging"]
    L = [f"# Stage 1 summary: {run['label']}", ""]
    L.append(f"- Target: **{run['target']}**   Line: **{run['line']}**   "
             f"Continuum subtraction: **{'none (nocontsub)' if run['nocontsub'] else 'uvcontsub applied'}**")
    if flag and flag["observation"]:
        L.append(f"- Observation (as logged in flag summary): {flag['observation']}")
    if flag and flag["ms"]:
        L.append(f"- MS at time of job: `{flag['ms']}`")
    L.append(f"- Products dir: `{run['dir']}`")
    L.append("")
    L.append("## Attention")
    L.append("")
    L += [f"- {a}" for a in run["attention"]] or ["- Nothing flagged automatically."]
    L.append("")

    # ---- flagging
    L += ["## Flagging", ""]
    if flag is None:
        L += ["No flag summary found.", ""]
    else:
        L.append(md_table(["", "before %", "after %", "change %", "data points"],
                          [["Total", f"{flag['before']:.3f}", f"{flag['after']:.3f}",
                            f"{flag['change']:+.3f}", f"{flag['npts']:,}"]], ["l", "r", "r", "r", "r"]))
        L.append("")
        fc = run["flagcmds"]
        if fc["n_cmds"]:
            L.append(f"Hub-antenna flag commands applied: **{fc['n_cmds']}** "
                     f"({len(fc['scans'])} scans, {len(fc['antennas'])} antennas, "
                     f"SPW(s) {fmt_list(fc['spws'])}).")
            L.append(f"- Scans: {fmt_list(fc['scans'], 20)}")
            L.append(f"- Antennas: {fmt_list(fc['antennas'], 30)}")
            L.append(f"- Files: {', '.join(f'`{f}`' for f in run['flagcmd_files'])}")
        else:
            L.append("No hub-antenna flag commands were applied (no hub antennas found).")
        L.append("")
        for sec in ("SPW", "field"):
            tab = flag["sections"].get(sec)
            if tab:
                L += [f"### Per {sec}", "", md_table([sec, "before %", "after %", "delta %"],
                                                     [[k, f"{b:.3f}", f"{a:.3f}", f"{d:.3f}"]
                                                      for k, (b, a, d) in tab.items()],
                                                     ["l", "r", "r", "r"]), ""]
        ants = flag["sections"].get("antenna")
        if ants:
            changed = [(k, v) for k, v in ants.items() if v[2] > 0]
            worst = sorted(ants.items(), key=lambda kv: -kv[1][1])[:8]
            L += ["### Per antenna", "",
                  f"{len(ants)} antennas; {len(changed)} had flags added. Highest flagged after flagging:", "",
                  md_table(["antenna", "before %", "after %", "delta %"],
                           [[k, f"{b:.3f}", f"{a:.3f}", f"{d:.3f}"] for k, (b, a, d) in worst],
                           ["l", "r", "r", "r"]), ""]
            if changed:
                L += ["Antennas with flags added: " +
                      ", ".join(f"{k} (+{v[2]:.3f}%)" for k, v in
                                sorted(changed, key=lambda kv: -kv[1][2])), ""]
        scans = flag["sections"].get("scan")
        if scans:
            changed = [(k, v) for k, v in scans.items() if v[2] > 0]
            top = sorted(scans.items(), key=lambda kv: -kv[1][1])[:8]
            L += ["### Per scan", "",
                  f"{len(scans)} scans; {len(changed)} had flags added. Most flagged after flagging:", "",
                  md_table(["scan", "before %", "after %", "delta %"],
                           [[k, f"{b:.3f}", f"{a:.3f}", f"{d:.3f}"] for k, (b, a, d) in top],
                           ["l", "r", "r", "r"]), ""]
            if changed:
                L += ["Scans with flags added: " +
                      ", ".join(f"{k} (+{v[2]:.2f}%)" for k, v in
                                sorted(changed, key=lambda kv: -kv[1][2])[:25]) +
                      (" ..." if len(changed) > 25 else ""), ""]

    # ---- imaging
    L += ["## Test imaging", ""]
    if imaging is None:
        L += ["No test imaging summary found.", ""]
    else:
        L.append(f"SPW {imaging['spw']}, channel {imaging['chan']} of {imaging['nchan']} "
                 f"({imaging['freq_mhz']:.3f} MHz), natural weighting, datacolumn={imaging['datacolumn']}.")
        L.append(f"Images: **{run['n_ok']} ok**, {run['n_fail_expected']} failed with no data "
                 f"in the imaged SPW/config, **{run['n_fail_unexplained']} failed unexplained**.")
        L.append("")
        L.append(md_table(["config", "obsids", "p95 baseline (m)", "median resid rms (mJy/beam)"],
                          [[c, fmt_list(v["obsids"], 8), v["p95"],
                            f"{run['median_rms'][c] * 1e3:.2f}" if c in run["median_rms"] else "-"]
                           for c, v in imaging["configs"].items()], ["l", "l", "r", "r"]))
        L.append("")
        trows = []
        for r in imaging["rows"]:
            if r["ok"] and "rms" in r:
                trows.append([r["field"], r["config"], f"{r['cell']:.2f}", r["imsize"],
                              f"{r['psf_max']:.3f}", f"{r['psf_min']:.3f}", f"{r['rms'] * 1e3:.3f}",
                              f"{r['resid_peak_snr']:.1f}", "ok"])
            elif r["ok"]:
                trows.append([r["field"], r["config"], f"{r['cell']:.2f}", r["imsize"],
                              "-", "-", f"{r['resid_rms'] * 1e3:.3f}", "-", f"ok ({r['note']})"])
            else:
                trows.append([r["field"], r["config"], f"{r['cell']:.2f}", r["imsize"], "-", "-", "-", "-",
                              f"FAILED, {r['note']}"])
        L.append(md_table(["field", "cfg", "cell (\")", "imsize", "PSF max", "PSF min",
                           "resid rms (mJy/bm)", "resid peak/rms", "status"], trows,
                          ["l", "l", "r", "r", "r", "r", "r", "r", "l"]))
        L.append("")
        if gallery_name:
            L += ["### Gallery", "", f"![{run['label']} gallery]({gallery_name})", ""]
    with open(path, "w") as fh:
        fh.write("\n".join(L))


def overview_rows(runs):
    rows = []
    for r in runs:
        f, im = r["flag"], r["imaging"]
        row = OrderedDict(
            run=r["label"], target=r["target"], line=r["line"],
            contsub="no" if r["nocontsub"] else "yes",
            flagged_before_pct=f["before"] if f else "", flagged_after_pct=f["after"] if f else "",
            flag_change_pct=f["change"] if f else "", n_flag_cmds=r["flagcmds"]["n_cmds"],
            imaged_spw=im["spw"] if im else "", imaged_freq_mhz=im["freq_mhz"] if im else "",
            n_images=len(im["rows"]) if im else 0, n_ok=r["n_ok"],
            n_fail_nodata=r["n_fail_expected"], n_fail_unexplained=r["n_fail_unexplained"],
        )
        for cfg in "ABCD":
            row[f"median_rms_{cfg}_mjy"] = (f"{r['median_rms'][cfg] * 1e3:.3f}" if cfg in r["median_rms"] else "")
        row["n_attention"] = len(r["attention"])
        rows.append(row)
    return rows


def write_overview(runs, outdir, data_dir):
    rows = overview_rows(runs)
    with open(os.path.join(outdir, "overview.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    L = ["# Stage 1 products overview", "", f"Data dir: `{data_dir}`   Runs: **{len(runs)}**", ""]
    tot_img = sum(r["n_ok"] + r["n_fail_expected"] + r["n_fail_unexplained"] for r in runs)
    L.append(f"Test images: {sum(r['n_ok'] for r in runs)} ok / "
             f"{sum(r['n_fail_expected'] for r in runs)} failed with no data in the imaged SPW or config / "
             f"{sum(r['n_fail_unexplained'] for r in runs)} failed unexplained (of {tot_img}).")
    L.append("")

    # runs where imaging is mostly empty
    bad = [r for r in runs if r["imaging"] and r["imaging"]["rows"]
           and r["n_ok"] < 0.5 * len(r["imaging"]["rows"])]
    if bad:
        causes = defaultdict(int)
        for r in bad:
            for row in r["imaging"]["rows"]:
                if not row["ok"]:
                    causes[row["note"]] += 1
        L += ["## Main finding", "",
              f"**{len(bad)} runs produced test images for fewer than half of the field/config combinations** "
              f"({', '.join(sorted({r['line'] for r in bad}))}). Imaged SPW(s): "
              f"{', '.join(sorted({r['imaging']['spw'] for r in bad}))}. Most common causes of the failed "
              "images (see each run summary for the per-image reason):", ""]
        L += [f"- {n} images: {c}" for c, n in sorted(causes.items(), key=lambda kv: -kv[1])[:3]]
        L += ["", "The affected runs:", ""]
        L += [f"- `{r['label']}`: {r['n_ok']}/{len(r['imaging']['rows'])} images ok" for r in bad]
        L.append("")

    L += ["## Runs", ""]
    hdr = ["run", "flagged before %", "after %", "added %", "flag cmds", "images ok/total",
           "failed (no data/unexpl.)", "median rms by config (mJy/bm)", "attention"]
    trows = []
    for r, o in zip(runs, rows):
        trows.append([f"[{r['label']}]({r['label']}_summary.md)",
                      f"{o['flagged_before_pct']}", f"{o['flagged_after_pct']}", f"{o['flag_change_pct']}",
                      o["n_flag_cmds"], f"{o['n_ok']}/{o['n_images']}",
                      f"{o['n_fail_nodata']}/{o['n_fail_unexplained']}",
                      " ".join(f"{c}:{o[f'median_rms_{c}_mjy']}" for c in "ABCD" if o[f"median_rms_{c}_mjy"]) or "-",
                      o["n_attention"]])
    L.append(md_table(hdr, trows, ["l"] + ["r"] * 8))
    L.append("")

    L += ["## Attention items by run", ""]
    for r in runs:
        if r["attention"]:
            L.append(f"### {r['label']}")
            L += [f"- {a}" for a in r["attention"]]
            L.append("")
    with open(os.path.join(outdir, "overview.md"), "w") as fh:
        fh.write("\n".join(L))


# =========================================================================================
# main
# =========================================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("data_dir", nargs="?", default=DEFAULT_DATA_DIR)
    ap.add_argument("-o", "--outdir", default=None)
    ap.add_argument("--match", default="*", help="glob on the run label, e.g. 'ic1613_*oh*'")
    ap.add_argument("--no-gallery", action="store_true")
    args = ap.parse_args()

    outdir = args.outdir or os.path.join(args.data_dir, "stage1_qa_reports")
    os.makedirs(outdir, exist_ok=True)
    dirs = sorted(d for d in glob.glob(os.path.join(args.data_dir, "*_products")) if os.path.isdir(d)
                  and fnmatch.fnmatch(run_label(d), args.match))
    if not dirs:
        sys.exit(f"No *_products directories matching '{args.match}' in {args.data_dir}")

    runs = []
    for d in dirs:
        print(f"[{len(runs) + 1}/{len(dirs)}] {run_label(d)}", flush=True)
        run = analyse_run(d)
        gallery_name = None
        if not args.no_gallery:
            gname = f"{run['label']}_gallery.png"
            try:
                if make_gallery(run, os.path.join(outdir, gname)):
                    gallery_name = gname
            except Exception as exc:
                print(f"  gallery failed: {exc}")
                run["attention"].append(f"gallery could not be made: {exc}")
        write_run_report(run, os.path.join(outdir, f"{run['label']}_summary.md"), gallery_name)
        runs.append(run)
    write_overview(runs, outdir, args.data_dir)
    print(f"Wrote {len(runs)} run reports + overview.md/csv to {outdir}")


if __name__ == "__main__":
    main()
