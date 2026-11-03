"""Tar helpers for the final staging pipeline (standard library only).

Used by job1_prep_and_test.sh (extract) and job2_archive_and_cleanup.sh (tar).
"""

import os
import tarfile


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
