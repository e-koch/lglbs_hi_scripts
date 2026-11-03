"""Shared helpers for the final staging pipeline (CASA-agnostic, except where noted).

Imported by run_uvcontsub.py / run_test_imaging.py (run under `casa -c`, which
defines __file__ so callers can sys.path.insert their own directory) and usable
from plain python.
"""

import os
import sys
import tarfile

import numpy as np

STAGED_DIR = os.environ.get(
    "STAGED_DIR",
    "/home/ekoch/projects/rrg-eros-ab/ekoch/VLAXL/staged_measurement_sets")
SCRATCH_ROOT = os.environ.get("IMAGING_SCRATCH", "/home/ekoch/scratch/VLAXL_imaging")


def script_args(script_name):
    """Arguments following `script_name` on a `casa -c script args` command line."""
    idx = max(i for i, a in enumerate(sys.argv) if a.endswith(script_name))
    return sys.argv[idx + 1:]


# ---------------------------------------------------------------------------
# Tar helpers
# ---------------------------------------------------------------------------

def extract_ms_tar(tar_path, dest_dir):
    """Extract an MS tar into dest_dir, stripping any leading path before '*.ms'.

    Staged tars store absolute-ish paths (home/ekoch/scratch/.../name.ms/...).
    Returns the path of the extracted MS directory.
    """
    ms_name = None
    with tarfile.open(tar_path, "r") as tf:
        for member in tf:
            parts = member.name.strip("/").split("/")
            ms_idx = next((i for i, p in enumerate(parts) if p.endswith(".ms")), None)
            if ms_idx is None:
                continue
            ms_name = parts[ms_idx]
            member.name = "/".join(parts[ms_idx:])
            tf.extract(member, path=dest_dir)
    if ms_name is None:
        raise RuntimeError(f"No *.ms directory found in {tar_path}")
    return os.path.join(dest_dir, ms_name)


def make_ms_tar(ms_path, tar_path):
    """Create an uncompressed tar of a directory MS with a relative '<name>.ms/' prefix."""
    ms_path = ms_path.rstrip("/")
    with tarfile.open(tar_path, "w") as tf:
        tf.add(ms_path, arcname=os.path.basename(ms_path), recursive=True)


def verify_tar(tar_path, ms_name):
    """Raise if tar can't be listed fully or lacks the MS table."""
    n = 0
    has_table = False
    with tarfile.open(tar_path, "r") as tf:
        for m in tf:
            n += 1
            if m.name == f"{ms_name}/table.dat":
                has_table = True
    if not has_table:
        raise RuntimeError(f"{tar_path}: {ms_name}/table.dat not found ({n} members)")
    return n


# ---------------------------------------------------------------------------
# VLA configuration identification (absolute, not rank-based)
# ---------------------------------------------------------------------------

# Nominal maximum baselines (m). The 95th-percentile uv distance of a track is
# ~0.7 x the maximum baseline, so compare p95 against 0.7 x nominal.
NOMINAL_MAX_BASELINE = {"A": 36.4e3, "B": 11.1e3, "C": 3.4e3, "D": 1.03e3}
P95_FRACTION = 0.7


def config_from_p95(p95):
    """Map a p95 uv distance (m) to A/B/C/D by the log-midpoint between configs."""
    order = ["A", "B", "C", "D"]
    nominal = {k: P95_FRACTION * NOMINAL_MAX_BASELINE[k] for k in order}
    for hi, lo in zip(order[:-1], order[1:]):
        if p95 >= np.sqrt(nominal[hi] * nominal[lo]):
            return hi
    return "D"


def obsid_p95_baselines(tb, vis, chunk_size=500000, uvmin=10.0, uvmax=1.0e5, nbins=400):
    """Per-OBSERVATION_ID p95 uv distance (m) from unflagged rows (FLAG_ROW only).

    `tb` is the CASA table tool. Uses a log-binned histogram so memory is bounded.
    """
    edges = np.logspace(np.log10(uvmin), np.log10(uvmax), nbins + 1)
    hist = {}
    tb.open(vis)
    nrows = tb.nrows()
    for start in range(0, nrows, chunk_size):
        n = min(chunk_size, nrows - start)
        obsid = tb.getcol("OBSERVATION_ID", startrow=start, nrow=n)
        uvw = tb.getcol("UVW", startrow=start, nrow=n)
        frow = tb.getcol("FLAG_ROW", startrow=start, nrow=n)
        uvd = np.hypot(uvw[0], uvw[1])
        good = (~frow) & (uvd > 0)
        idx = np.clip(np.digitize(uvd, edges) - 1, 0, nbins - 1)
        for oid in np.unique(obsid[good]):
            sel = good & (obsid == oid)
            h = hist.setdefault(int(oid), np.zeros(nbins, dtype=np.int64))
            h += np.bincount(idx[sel], minlength=nbins)
    tb.close()

    out = {}
    for oid, h in hist.items():
        cum = np.cumsum(h) / h.sum()
        i = min(int(np.searchsorted(cum, 0.95)), nbins - 1)
        lo = cum[i - 1] if i > 0 else 0.0
        frac = 0.0 if cum[i] == lo else (0.95 - lo) / (cum[i] - lo)
        out[oid] = edges[i] + frac * (edges[i + 1] - edges[i])
    return out


def identify_configs(tb, vis):
    """Return ({config: [obsids]}, {config: p95 baseline (m)}, {obsid: p95})."""
    p95 = obsid_p95_baselines(tb, vis)
    groups = {}
    for oid, v in sorted(p95.items()):
        groups.setdefault(config_from_p95(v), []).append(oid)
    group_p95 = {c: float(np.max([p95[o] for o in oids])) for c, oids in groups.items()}
    return groups, group_p95, p95


# ---------------------------------------------------------------------------
# Imaging-size helper
# ---------------------------------------------------------------------------

def good_imsize(n, max_size=4096):
    """Smallest even number >= n whose only prime factors are 2, 3, 5 (capped)."""
    n = int(min(max(n, 64), max_size))
    while n <= max_size:
        m = n
        for f in (2, 3, 5):
            while m % f == 0:
                m //= f
        if m == 1 and n % 2 == 0:
            return n
        n += 1
    return max_size
