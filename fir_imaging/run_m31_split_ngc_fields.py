'''
Split the NGC185 / NGC205 fields out of the staged per-config M31 MSs.

The 14A-235 D-config tracks also observed NGC185 and NGC205. The early per-config M31 staging
selected every field (ms_file_key: `all`), so those fields ended up in m31_D_<line>[_nocontsub].ms.tar
and from there in the C+D and A+B+C+D concats. The staging key now selects `M31*` only, so newly
staged tars are clean; this script repairs the tars staged before that fix.

For every staged tar  <staged>/m31_<cfg>_<line>[_nocontsub].ms.tar  (cfg from M31_SPLIT_CONFIGS,
default D; the concatenated m31_C+D_* / m31_A+B+C+D_* tars never match) this script:
1. Extracts the MS to local scratch (one tar at a time).
2. Lists the fields that have data. NGC185/NGC205 absent -> already clean, tar left untouched.
   A field that is neither M31* nor NGC185/NGC205 -> error, tar left untouched.
3. split()s the M31* fields and each NGC field into their own MS (datacolumn='all',
   keepflags=True: nothing is dropped, flags are kept exactly).
4. Checks that the row counts add up to the original and that no MAIN column was lost.
5. Writes ngc185_<cfg>_<line>[_nocontsub].ms.tar and ngc205_<...>.ms.tar into the staged dir,
   then replaces the m31 tar with the M31*-only tar (each written as <name>.tar.tmp and renamed
   when complete, and the m31 tar is replaced LAST so an interrupted run can simply be repeated).
6. Cleans up the local files.

Safe to re-run: a tar that no longer contains NGC fields is skipped.

Environment variables (set by job_split_m31_ngc_fields.sh, which passes its own optional env):
    M31_SPLIT_LOCAL_DIR  per-job scratch dir (required to be set by the job script)
    M31_SPLIT_CONFIGS    comma list of configs to process (default D)
    M31_SPLIT_LINES      comma list of line names, exact match (default: all lines found)
    M31_SPLIT_STAGED_DIR staged_measurement_sets dir (default: the project-space one; point it at a
                         test copy to try this out)
    M31_SPLIT_DRY_RUN    1 = only read the FIELD table of each tar and report which NGC field names
                         it lists; nothing is extracted/written. (The FIELD table can list fields that
                         have no data; the real run checks for scans and skips those tars.)
'''

import glob
import os
import shutil
import sys
import tarfile

from casatasks import split  # type: ignore  # only available inside CASA
from casatools import msmetadata, table  # type: ignore  # only available inside CASA

# pipeline_common.py (tar helpers) lives in final_staging_pipeline/
pipe_dir = os.environ.get('M31_SPLIT_PIPE_DIR',
                          '/home/ekoch/lglbs_hi_scripts/fir_imaging/final_staging_pipeline')
sys.path.insert(0, pipe_dir)
from pipeline_common import extract_ms_tar, make_ms_tar, verify_tar  # noqa: E402

tb = table()
msmd = msmetadata()

targ = 'm31'
ngc_targets = ['ngc185', 'ngc205']     # field-name prefix (case-insensitive) -> output target name

staged_dir = os.environ.get('M31_SPLIT_STAGED_DIR',
                            '/home/ekoch/projects/rrg-eros-ab/ekoch/VLAXL/staged_measurement_sets/')
local_dir = os.environ.get('M31_SPLIT_LOCAL_DIR', '.')
os.makedirs(local_dir, exist_ok=True)


def _env_list(name, default=None):
    val = os.environ.get(name, '').strip()
    return [v.strip() for v in val.split(',') if v.strip()] if val else default


configs = _env_list('M31_SPLIT_CONFIGS', ['D'])
lines = _env_list('M31_SPLIT_LINES')
dry_run = os.environ.get('M31_SPLIT_DRY_RUN', '0').strip().lower() in ('1', 'true', 'yes', 'y')

print(f"Staged dir: {staged_dir}")
print(f"Local working dir: {os.path.abspath(local_dir)}")
print(f"configs={configs} lines={lines or 'all'} dry_run={dry_run}")


def classify(name):
    """'m31', one of ngc_targets, or None for a field we do not know how to place."""
    up = name.upper()
    if up.startswith('M31'):
        return targ
    for t in ngc_targets:
        if up.startswith(t.upper()):
            return t
    return None


def line_of(tar_name, cfg):
    """Line (+ '_nocontsub' suffix kept separately) from m31_<cfg>_<line>[_nocontsub].ms.tar."""
    stem = tar_name[:-len('.ms.tar')][len(f"{targ}_{cfg}_"):]
    nocs = stem.endswith('_nocontsub')
    return (stem[:-len('_nocontsub')] if nocs else stem), nocs


def field_names_from_tar(tar_path, dest):
    """Read FIELD/NAME from a tar without extracting the main table (dry run)."""
    ms_dir = None
    with tarfile.open(tar_path, 'r') as tf:
        for member in tf:
            parts = member.name.replace('\\', '/').split('/')
            ms_idx = next((i for i, p in enumerate(parts) if p.endswith('.ms')), None)
            if ms_idx is None or len(parts) <= ms_idx + 1 or parts[ms_idx + 1] != 'FIELD':
                continue
            ms_dir = parts[ms_idx]
            member.name = '/'.join(parts[ms_idx:])
            tf.extract(member, path=dest)
    if ms_dir is None:
        raise RuntimeError(f"{tar_path}: no FIELD subtable found")
    tb.open(os.path.join(dest, ms_dir, 'FIELD'))
    names = list(tb.getcol('NAME'))
    tb.close()
    shutil.rmtree(os.path.join(dest, ms_dir))
    return names


def write_tar(ms_path, dest_tar):
    """Tar ms_path (whose dir name must be <tar name minus .tar>) to dest_tar via .tmp + rename."""
    tmp = dest_tar + '.tmp'
    if os.path.exists(tmp):
        os.remove(tmp)
    try:
        make_ms_tar(ms_path, tmp)
        verify_tar(tmp, os.path.basename(ms_path.rstrip('/')))
        os.replace(tmp, dest_tar)
    except BaseException:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise


def process_tar(tar_path, cfg):
    """Returns 'clean', 'split' or 'dry'. Raises on anything unexpected (tar left untouched)."""
    tar_name = os.path.basename(tar_path)
    base = tar_name[:-len('.tar')]                      # m31_D_hilores_nocontsub.ms

    if dry_run:
        names = field_names_from_tar(tar_path, local_dir)
        ngc = sorted({n for n in names if classify(n) in ngc_targets})
        other = sorted({n for n in names if classify(n) is None})
        print(f"  FIELD table: {len(names)} fields; NGC fields listed: {ngc or 'none'}"
              + (f"; UNRECOGNISED: {other}" if other else ''))
        return 'dry'

    print("  Extracting")
    ms = extract_ms_tar(tar_path, local_dir)
    assert os.path.basename(ms) == base, f"MS dir {ms} does not match tar name {tar_name}"

    try:
        msmd.open(ms)
        try:
            names = list(msmd.fieldnames())
            with_data = [i for i, n in enumerate(names) if len(msmd.scansforfield(i)) > 0]
        finally:
            msmd.close()

        groups = {}
        for i in with_data:
            tgt = classify(names[i])
            if tgt is None:
                raise RuntimeError(f"field {i} '{names[i]}' is neither M31* nor "
                                   f"{'/'.join(n.upper() for n in ngc_targets)}; not touching this tar")
            groups.setdefault(tgt, []).append(i)
        print("  Fields with data: " + ', '.join(f"{t}: {len(v)}" for t, v in groups.items()))

        if not any(t in groups for t in ngc_targets):
            print("  No NGC fields: already clean, leaving the tar untouched")
            return 'clean'
        if targ not in groups:
            raise RuntimeError("no M31* fields found; not touching this tar")

        tb.open(ms)
        n_orig = tb.nrows()
        cols_orig = set(tb.colnames())
        tb.close()

        # ---- split every group into its own MS: <target>_<cfg>_<line>[_nocontsub].ms
        out = {}
        for tgt, ids in groups.items():
            out_ms = os.path.join(local_dir, base.replace(f"{targ}_", f"{tgt}_", 1))
            if os.path.exists(out_ms):
                shutil.rmtree(out_ms)
            print(f"  split -> {os.path.basename(out_ms)} ({len(ids)} fields)")
            split(vis=ms, outputvis=out_ms, field=','.join(str(i) for i in ids),
                  datacolumn='all', keepflags=True)
            out[tgt] = out_ms

        # ---- checks: nothing lost
        n_new = 0
        for tgt, out_ms in out.items():
            tb.open(out_ms)
            n = tb.nrows()
            missing = cols_orig - set(tb.colnames())
            tb.close()
            print(f"    {os.path.basename(out_ms)}: {n} rows")
            if missing:
                raise RuntimeError(f"{os.path.basename(out_ms)} lost MAIN columns {sorted(missing)}")
            n_new += n
        if n_new != n_orig:
            raise RuntimeError(f"row count mismatch: original {n_orig}, split pieces sum to {n_new}")
        print(f"  Row counts OK: {n_orig} = {n_new}")

        # ---- write: NGC tars first, the (cleaned) m31 tar last
        for tgt in [t for t in out if t != targ] + [targ]:
            dest = os.path.join(staged_dir, os.path.basename(out[tgt]) + '.tar')
            print(f"  Writing {dest}")
            write_tar(out[tgt], dest)
        return 'split'

    finally:
        for path in [ms] + [os.path.join(local_dir, base.replace(f"{targ}_", f"{t}_", 1))
                            for t in [targ] + ngc_targets]:
            if os.path.exists(path):
                shutil.rmtree(path)


results = {}
failed = []
for cfg in configs:
    tars = sorted(glob.glob(os.path.join(staged_dir, f"{targ}_{cfg}_*.ms.tar")))
    # m31_D_* does not match the concatenated m31_C+D_* / m31_A+B+C+D_* tars
    if lines:
        tars = [t for t in tars if line_of(os.path.basename(t), cfg)[0] in lines]
    print(f"\n##### config {cfg}: {len(tars)} tar(s)")
    for tar_path in tars:
        print(f"=== {os.path.basename(tar_path)} ===")
        try:
            results[tar_path] = process_tar(tar_path, cfg)
        except Exception as exc:            # keep going; report at the end
            print(f"  ERROR: {exc}")
            results[tar_path] = 'ERROR'
            failed.append(tar_path)

print("\n##### Summary")
for tar_path, res in results.items():
    print(f"  {res:>6}  {os.path.basename(tar_path)}")
if failed:
    raise SystemExit(f"{len(failed)} tar(s) failed; see the log above")
