
'''
The M31 staging was large enough that it was easier to do per config.

Now we want to put it back together.

For each line product, this script:
1. Copies each per-config tar from staged_measurement_sets to local scratch.
2. Untars each, removes the local tar.
3. CASA-concats the configs into a combined MS: A+B+C+D or C+D, chosen with the
   M31_CONCAT_COMBO env var (ABCD or CD; set by job_submit_m31_concat.sh <combo>).
   Work happens in M31_CONCAT_LOCAL_DIR (a per-job scratch dir).
4. Tars the combined MS using the standard naming convention, written directly into
   staged_measurement_sets (as <name>.tar.tmp, renamed to <name>.tar only when complete, so an
   interrupted job never leaves a partial tar under the real name; no local copy of the tar).
5. Cleans up all local files for that product.

Overwrite behaviour (overwrite=True by default; set M31_CONCAT_OVERWRITE=0 to change):
- overwrite=True : an existing <name>.tar in staged_measurement_sets is removed just before the
                   new tar is written (after the concat succeeded).
- overwrite=False: a line whose combined tar already exists is skipped entirely (no copy,
                   extract or concat).
'''

import os
import shutil
import tarfile

import numpy as np
from casatasks import concat  # type: ignore  # only available inside CASA
from casatools import table  # type: ignore  # only available inside CASA

tb = table()
freqtol_chanfrac = 0.2  # phangs_imaging_pipeline default

targ = 'm31'

# Which configs to combine, set by job_submit_m31_concat.sh via M31_CONCAT_COMBO.
# ABCD -> m31_A+B+C+D_<line>.ms ; CD -> m31_C+D_<line>.ms
combo_options = {'ABCD': ['A', 'B', 'C', 'D'],
                 'CD': ['C', 'D']}
combo = os.environ.get('M31_CONCAT_COMBO', 'ABCD').upper().replace('+', '')
if combo not in combo_options:
    raise ValueError(f"M31_CONCAT_COMBO={combo!r} not in {list(combo_options)}")
configs = combo_options[combo]
combo_label = '+'.join(configs)
print(f"Combining configs {combo_label}")

all_line_products = ['hilores', 'himidres', 'hi21cm_0p8kms', 'hi21cm',
                     'oh1612', 'oh1720', 'oh1665', 'oh1667']

data_path = "/home/ekoch/scratch/VLAXL_imaging/"
export_path = "/home/ekoch/projects/rrg-eros-ab/ekoch/VLAXL/staged_measurement_sets/"
# Per-job scratch dir (set by the submit script) so parallel jobs don't collide.
local_dir = os.environ.get('M31_CONCAT_LOCAL_DIR', '.')
os.makedirs(local_dir, exist_ok=True)
print(f"Local working dir: {os.path.abspath(local_dir)}")

# Replace an existing combined tar (True) or skip that line (False). Set by the submit command via
# M31_CONCAT_OVERWRITE (e.g. `M31_CONCAT_OVERWRITE=0 sbatch job_submit_m31_concat.sh CD`).
overwrite = os.environ.get('M31_CONCAT_OVERWRITE', '1').strip().lower() not in ('0', 'false', 'no', 'n')
print(f"overwrite={overwrite}")

for this_line in all_line_products:
    print(f"=== Processing {this_line} ===")

    combined_tar_name = f"{targ}_{combo_label}_{this_line}.ms.tar"
    dest_tar = os.path.join(export_path, combined_tar_name)
    tmp_tar = dest_tar + ".tmp"
    if os.path.exists(tmp_tar):
        # Left behind by an interrupted job; never a complete tar.
        print(f"  Removing stale partial tar {tmp_tar}")
        os.remove(tmp_tar)
    if os.path.exists(dest_tar) and not overwrite:
        print(f"  {combined_tar_name} already exists and overwrite=False; skipping {this_line}\n")
        continue

    # Steps 1+2: Copy each per-config tar, untar, remove local tar
    per_config_ms = []
    for cfg in configs:
        ms_name = f"{targ}_{cfg}_{this_line}.ms"
        tar_name = ms_name + ".tar"
        src = os.path.join(export_path, tar_name)
        dst = os.path.join(local_dir, tar_name)

        print(f"  Copying {tar_name}")
        shutil.copy2(src, dst)

        print(f"  Extracting {tar_name}")
        # The tar contains absolute paths (e.g. home/ekoch/scratch/.../m31_A_hilores.ms/...).
        # Strip all leading components up to and including the .ms directory so
        # the MS lands directly in local_dir.
        with tarfile.open(dst, "r") as tfile:
            for member in tfile.getmembers():
                parts = member.name.replace("\\", "/").split("/")
                # Find the component that ends with .ms and drop everything before it
                ms_idx = next((i for i, p in enumerate(parts) if p.endswith(".ms")), None)
                if ms_idx is None:
                    continue
                member.name = "/".join(parts[ms_idx:])
                tfile.extract(member, path=local_dir)
        os.remove(dst)

        per_config_ms.append(os.path.join(local_dir, ms_name))

    # Step 3: CASA concat all 4 configs
    combined_ms_name = f"{targ}_{combo_label}_{this_line}.ms"
    combined_ms = os.path.join(local_dir, combined_ms_name)
    print(f"  Concatenating into {combined_ms_name}")
    # Tolerance for merging SPWs: same rule as phangs_imaging_pipeline
    # (casaVisRoutines.py: freqtol = freqtol_chanfrac * final chan width,
    # freqtol_chanfrac=0.2). Tiny chan-width / frequency offsets between configs
    # then merge into one SPW instead of creating a new SPW per config.
    max_chan_hz = 0.
    for ms in per_config_ms:
        tb.open(os.path.join(ms, 'SPECTRAL_WINDOW'))
        max_chan_hz = max(max_chan_hz, float(np.abs(tb.getcol('CHAN_WIDTH')).max()))
        tb.close()
    freqtol = f"{freqtol_chanfrac * max_chan_hz}Hz"
    print(f"  Using freqtol={freqtol} ({freqtol_chanfrac} x max chan width {max_chan_hz} Hz)")
    concat(vis=per_config_ms, concatvis=combined_ms, freqtol=freqtol,
           dirtol='10arcsec',  # tiny fraction of the L-band primary beam
           copypointing=False)

    # Step 4: Tar the combined MS directly into the project space (no local copy of the tar).
    # Written as <name>.tar.tmp (not matched by the *.ms.tar globs) and renamed once complete.
    if os.path.exists(dest_tar):   # only reachable with overwrite=True
        print(f"  Removing existing {combined_tar_name} (overwrite=True)")
        os.remove(dest_tar)
    print(f"  Creating tar {dest_tar}")
    try:
        with tarfile.open(tmp_tar, "w") as tfile:
            tfile.add(combined_ms, recursive=True)
        os.rename(tmp_tar, dest_tar)
    except BaseException:
        if os.path.exists(tmp_tar):
            os.remove(tmp_tar)
        raise

    # Step 5: Clean up all local MS directories for this product
    print(f"  Cleaning up local files for {this_line}")
    for cfg in configs:
        ms_path = os.path.join(local_dir, f"{targ}_{cfg}_{this_line}.ms")
        if os.path.exists(ms_path):
            shutil.rmtree(ms_path)
    if os.path.exists(combined_ms):
        shutil.rmtree(combined_ms)

    print(f"=== Done with {this_line} ===\n")
