
'''
Stage the wide-field C+D config data for the dwarfs.

Produces staging for:
- 4 OH lines
- HI at 4.2 km/s (10x native)
- HI at 2.1 km/s (5x native)
- HI at native 0.4 km/s
- Foreground HI at native 0.4 km/s (where appropriate)


Run in CASA from the command line.
Command line argument is the galaxy name.

'''

##############################################################################
# Load routines, initialize handlers
##############################################################################

import os
import sys
from importlib import reload

# Locate the master key
#key_file = '/data/tycho/0/leroy.42/reduction/vla/lglbs_pipeline_configs/lglbs_keys/master_key.txt'
key_file = "/home/erickoch/lglbs_hi_scripts/lglbs_keys/master_key.txt"

data_path = "/reduction10/erickoch/LGLBS/line_imaging/"
export_path = f"{data_path}/staged_measurement_sets/"

# Set the logging
from phangsPipeline import phangsLogger as pl
# reload(pl)
pl.setup_logger(level='DEBUG', logfile=None)

from phangsPipeline import handlerKeys as kh
from phangsPipeline import handlerVis as uvh
# from phangsPipeline import handlerImaging as imh
# from phangsPipeline import handlerPostprocess as pph

# Initialize key handlers
this_kh = kh.KeyHandler(master_key = key_file)
this_uvh = uvh.VisHandler(key_handler = this_kh)

# Make missing directories
this_kh.make_missing_directories(imaging=True,derived=True,postprocess=True,release=True)

##############################################################################
# Set up what we do this run
##############################################################################

this_targ = 'm33'


# all_configs = ['C+D', 'A+B+C+D', 'B+C+D', 'A+B+C']
# all_configs = ['C+D'] #, 'B+C+D']
all_configs = ['D']

this_uvh.set_targets(only=[this_targ])
this_uvh.set_interf_configs(only=all_configs)

# all_line_products = ['hi21cm', 'hi21cm_0p8kms', 'hilores', 'himidres',
#                      'oh1612', 'oh1720', 'oh1665', 'oh1667']

# all_line_products = ['oh1612', 'oh1720', 'oh1665', 'oh1667']
all_line_products = ['hi21cm_0p8kms', 'hilores', 'himidres']

# this_uvh.set_line_products(only=all_line_products)

this_uvh.set_no_cont_products(True)

##############################################################################
# Run staging
##############################################################################

for this_line in all_line_products:
    print(f"Working on {this_line}")

    this_uvh.set_line_products(only=[this_line])

    this_uvh.loop_stage_uvdata(do_copy=True, do_contsub=True,
                                do_extract_line=False, do_extract_cont=False,
                                do_remove_staging=False, overwrite=True, strict_config=False)

    this_uvh.loop_stage_uvdata(do_copy=False, do_contsub=False,
                                do_extract_line=True, do_extract_cont=False,
                                do_remove_staging=False, overwrite=True, strict_config=False)

    # this_uvh.loop_stage_uvdata(do_copy=False, do_contsub=False,
    #                             do_extract_line=False, do_extract_cont=True,
    #                             do_remove_staging=False, overwrite=True, strict_config=False)

    # NOTE: READD THIS AFTER WLM CONTSUB TESTS ARE DONE
    this_uvh.loop_stage_uvdata(do_copy=False, do_contsub=False,
                                do_extract_line=False, do_extract_cont=False,
                                do_remove_staging=True, overwrite=True, strict_config=False)


# Now we'll tar up the MS files

import tarfile

for this_config in all_configs:
    for this_line in all_line_products:

        final_msname = f"{this_targ}_{this_config}_{this_line}.ms"

        full_final_msname = f"{data_path}/imaging/{this_targ}/{final_msname}"

        if not os.path.exists(full_final_msname):
            print(f"Unable to find an existing ms file: {full_final_msname}. Skipping tar")
            continue

        # Create a new tar file. Remove the old version if it already exists.
        tar_msname = f"{export_path}/{final_msname}.tar"
        if os.path.exists(tar_msname):
            print(f"Found existing {tar_msname}. Deleting and creating new tar file.")
            os.remove(tar_msname)

        with tarfile.open(tar_msname, "w") as tfile:
            tfile.add(full_final_msname, recursive=True)
