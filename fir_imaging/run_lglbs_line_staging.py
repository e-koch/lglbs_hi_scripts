
'''
Stage LGLBS line data

Produces staging for:
- 4 OH lines
- HI at 4.2 km/s (10x native)
- HI at 2.1 km/s (5x native)
- HI at native 0.4 km/s
- HI at 0.84 km/s

'''

##############################################################################
# Load routines, initialize handlers
##############################################################################

import os
import sys
import tarfile


# Locate the master key
key_file = "/home/ekoch/lglbs_hi_scripts/lglbs_keys/master_key_fir.txt"

data_path = "/home/ekoch/scratch/VLAXL_imaging/"
export_path = "/home/ekoch/projects/rrg-eros-ab/ekoch/VLAXL/staged_measurement_sets/"

# Set the logging
from phangsPipeline import phangsLogger as pl
# reload(pl)
pl.setup_logger(level='DEBUG', logfile=None)

from phangsPipeline import handlerKeys as kh
from phangsPipeline import handlerVis as uvh


# Initialize key handlers
this_kh = kh.KeyHandler(master_key = key_file)
this_uvh = uvh.VisHandler(key_handler = this_kh)

# Make missing directories
this_kh.make_missing_directories(imaging=True,
                                 derived=True,
                                 postprocess=True,
                                 release=True)

##############################################################################
# Set up what we do this run
##############################################################################

# this_targ = 'm31'
this_targ = sys.argv[-2].lower()

this_config = sys.argv[-1]

# Varying spatial coverage in A, B vs C, D. Specify central pointing only for these
dwarfs_targs = ['wlm', 'ic1613', 'ic10']

if this_targ in dwarfs_targs and "A" in this_config:
    # Grab only the central pointing for the dwarfs
    this_targ += 'ctr'

# all_configs = ['C+D', 'A+B+C+D', 'B+C+D', 'A+B+C']
# all_configs = ['C+D', 'D'] #, 'B+C+D']
# all_configs = ['A+B+C+D']

all_configs = [this_config]

this_uvh.set_targets(only=[this_targ])
this_uvh.set_interf_configs(only=all_configs)

# all_line_products = ['hi21cm', 'hi21cm_0p8kms', 'hilores', 'himidres',
#                      'oh1612', 'oh1720', 'oh1665', 'oh1667']

all_line_products = []
# all_line_products += ['oh1612', 'oh1720', 'oh1665', 'oh1667']
all_line_products += ['himidres']
all_line_products += ['hilores']
all_line_products += ['hi21cm_0p8kms']
all_line_products += ['hi21cm']

# this_uvh.set_line_products(only=all_line_products)

this_uvh.set_no_cont_products(True)

##############################################################################
# Run staging
##############################################################################

# Make versions with and without contsub

# for do_contsub in [False, True]:
for do_contsub in [True]:

    if not do_contsub:
        contsub_str = "_nocontsub"
    else:
        contsub_str = ""

    for this_line in all_line_products:
        print(f"Working on {this_line}")

        this_uvh.set_line_products(only=[this_line])

        this_uvh.loop_stage_uvdata(do_copy=True, do_contsub=do_contsub,
                                    do_extract_line=False, do_extract_cont=False,
                                    do_remove_staging=False, overwrite=True,
                                    strict_config=False,)

        this_uvh.loop_stage_uvdata(do_copy=False, do_contsub=False,
                                    do_extract_line=True, do_extract_cont=False,
                                    do_remove_staging=False, overwrite=True,
                                    strict_config=False,)

        this_uvh.loop_stage_uvdata(do_copy=False, do_contsub=False,
                                    do_extract_line=False, do_extract_cont=False,
                                    do_remove_staging=True, overwrite=True,
                                    strict_config=False,)


        # Now we'll tar up the MS files _and_ rename files if they are not contsub'd

        for this_config in all_configs:

            final_msname = f"{this_targ}_{this_config}_{this_line}.ms"

            full_final_msname = f"{data_path}/imaging/{this_targ}/{final_msname}"

            if not os.path.exists(full_final_msname):
                print(f"Unable to find an existing ms file: {full_final_msname}. Skipping tar")
                continue

            if not do_contsub:
                # Rename the MS folder to avoid overwriting the contsub version
                final_msname_nocontsub = f"{this_targ}_{this_config}_{this_line}{contsub_str}.ms"

                if os.path.exists(final_msname_nocontsub):
                    # Remove the old version:
                    os.system(f"rm -r {final_msname_nocontsub}")

                os.rename(final_msname, final_msname_nocontsub)

                final_msname = final_msname_nocontsub
                full_final_msname = f"{data_path}/imaging/{this_targ}/{final_msname}"


            # Create a new tar file. Remove the old version if it already exists.
            tar_msname = f"{export_path}/{final_msname}.tar"
            if os.path.exists(tar_msname):
                print(f"Found existing {tar_msname}. Deleting and creating new tar file.")
                os.remove(tar_msname)

            with tarfile.open(tar_msname, "w") as tfile:
                tfile.add(full_final_msname, recursive=True)
