"""Helpers for the final staging pipeline (standard library only).

Tar helpers are used by job1_prep_and_test.sh (extract) and job2_archive_and_cleanup.sh (tar).
The channel-selection helpers are used by run_test_imaging.py and find_weight_outliers.py.
"""

import csv
import os
import tarfile

# Fraction of each SPW's channels at each end that is treated as line-free. Same value as
# edge_frac in run_uvcontsub.py (the uvcontsub fit range).
EDGE_FRAC = 0.05


def edge_channel_count(nchan, frac=EDGE_FRAC):
    """Number of channels at each end of an SPW (same rounding as run_uvcontsub.py)."""
    return max(1, int(frac * nchan + 0.5))


def read_flag_vs_channel_csv(csv_path, nchan):
    """Flagged % per channel from a plot_flag_vs_channel.py CSV, or None if missing/unusable."""
    if not os.path.isfile(csv_path):
        return None
    try:
        with open(csv_path) as fh:
            pct = [float(row["flagged_pct"]) for row in csv.DictReader(fh)]
    except (OSError, ValueError, KeyError):
        return None
    return pct if len(pct) == nchan else None


def least_flagged_edge_channel(pct, frac=EDGE_FRAC, tol=0.01):
    """Pick the least-flagged channel within the outer `frac` of channels at each end.

    pct: flagged % per channel. Channels within `tol` percentage points of the minimum count as
    tied; the tie is broken in favour of the channel furthest from the band edge (edge channels
    are the least reliable after regridding). Returns (channel, flagged_pct).
    """
    n = len(pct)
    k = edge_channel_count(n, frac)
    cands = list(range(k)) + list(range(n - k, n))
    best = min(pct[c] for c in cands)
    near = [c for c in cands if pct[c] <= best + tol]
    chan = max(near, key=lambda c: (min(c, n - 1 - c), -c))
    return chan, pct[chan]


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
