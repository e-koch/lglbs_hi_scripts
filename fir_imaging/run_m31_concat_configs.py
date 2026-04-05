
'''
The M31 staging was large enough that it was easier to do per config.

Now we want to put it back together.

For each line product, this script:
1. Copies each per-config tar from staged_measurement_sets to local scratch.
2. Untars each, removes the local tar.
3. CASA-concats all 4 configs into a combined A+B+C+D MS.
4. Tars the combined MS using the standard naming convention.
5. Moves the combined tar back to staged_measurement_sets.
6. Cleans up all local files for that product.
'''

import os
import shutil
import subprocess
import tarfile

from casatasks import concat  # type: ignore  # only available inside CASA

targ = 'm31'
configs = ['A', 'B', 'C', 'D']
all_line_products = ['hilores', 'himidres', 'hi21cm_0p8kms', 'hi21cm',
                     'oh1612', 'oh1720', 'oh1665', 'oh1667']

data_path = "/home/ekoch/scratch/VLAXL_imaging/"
export_path = "/home/ekoch/projects/rrg-eros-ab/ekoch/VLAXL/staged_measurement_sets/"
local_dir = '.'
# local_dir = os.path.join(data_path, "imaging", targ)
# os.makedirs(local_dir, exist_ok=True)

for this_line in all_line_products:
    print(f"=== Processing {this_line} ===")

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
        # Tars store just the bare MS name, so extract into local_dir.
        subprocess.run(["tar", "-xf", dst, "-C", local_dir], check=True)
        os.remove(dst)

        per_config_ms.append(os.path.join(local_dir, ms_name))

    # Step 3: CASA concat all 4 configs
    combined_ms_name = f"{targ}_A+B+C+D_{this_line}.ms"
    combined_ms = os.path.join(local_dir, combined_ms_name)
    print(f"  Concatenating into {combined_ms_name}")
    concat(vis=per_config_ms, concatvis=combined_ms)

    # Step 4: Tar the combined MS
    combined_tar_name = combined_ms_name + ".tar"
    combined_tar = os.path.join(local_dir, combined_tar_name)
    print(f"  Creating tar {combined_tar_name}")
    with tarfile.open(combined_tar, "w") as tfile:
        tfile.add(combined_ms, recursive=True)

    # Step 5: Move combined tar to export_path
    dest_tar = os.path.join(export_path, combined_tar_name)
    print(f"  Moving tar to {export_path}")
    shutil.move(combined_tar, dest_tar)

    # Step 6: Clean up all local MS directories for this product
    print(f"  Cleaning up local files for {this_line}")
    for cfg in configs:
        ms_path = os.path.join(local_dir, f"{targ}_{cfg}_{this_line}.ms")
        if os.path.exists(ms_path):
            shutil.rmtree(ms_path)
    if os.path.exists(combined_ms):
        shutil.rmtree(combined_ms)

    print(f"=== Done with {this_line} ===\n")
