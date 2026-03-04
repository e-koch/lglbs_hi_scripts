
'''
Stage the MW foreground wide-field C+D config data for the dwarfs.
Note NGC6822 overlaps with the foreground so isn't processed separately.

This script is for the other 3 dwarfs.

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
key_file = "/home/erickoch/lglbs_hi_scripts/lglbs_mwforeground_keys/master_key.txt"


# Set directory for the pipeline and change to this directory
#pipedir = '/data/tycho/0/leroy.42/reduction/alma/phangs_imaging_scripts/'
#os.chdir(pipedir)

# Make sure we are inside CASA (modify this to use the command line version)
#sys.path.append(os.getcwd())
#casa_enabled = (sys.argv[0].endswith('start_casa.py'))
#if not casa_enabled:
#    print('Please run this script inside CASA!')
#    sys.exit()

# Set the logging
from phangsPipeline import phangsLogger as pl
reload(pl)
pl.setup_logger(level='DEBUG', logfile=None)

# Imports

#sys.path.insert(1, )
from phangsPipeline import handlerKeys as kh
from phangsPipeline import handlerVis as uvh
from phangsPipeline import handlerImaging as imh
from phangsPipeline import handlerPostprocess as pph

# Reloads for debugging
reload(kh)
reload(uvh)
reload(imh)
reload(pph)

# Initialize key handlers
this_kh = kh.KeyHandler(master_key = key_file)
this_uvh = uvh.VisHandler(key_handler = this_kh)

# Make missing directories
this_kh.make_missing_directories(imaging=True,derived=True,postprocess=True,release=True)

##############################################################################
# Set up what we do this run
##############################################################################

this_targ = sys.argv[-1]

all_targs = ['ic10', 'ic1613', 'wlm', 'ic10ctr', 'ic1613ctr', 'wlmctr']
all_targs = [f"{targ}-mwfore" for targ in all_targs]

if not this_targ in all_targs:
    raise ValueError(f"Cannot find target {this_targ} in list of target names: {all_targs}")


if 'ctr' in this_targ:
    # Only the center pointings
    print("This is a center pointing. Using all configs.")
    this_uvh.set_interf_configs(only=['B', 'A+B', 'B+C+D', 'A+B+C+D'])
else:
    # Full mosaic only in C, D configs
    print("This is a full mosaic. Using C and D configs.")
    this_uvh.set_interf_configs(only=['C+D', 'C', 'D'])

print(this_targ)

this_uvh.set_targets(only=[this_targ])
print(this_uvh.get_targets())

# Because we have hand-tuned the freq range to exclude per galaxy (so uvcontsub ignores
# both the MW and galaxy), the line products are keyed per target:
this_gal = this_targ.split("-mwfore")[0]
all_line_products = [f'hi21cm_{this_gal}',
                     f'hi21cm_1p2kms_{this_gal}',
                     f'hilores_{this_gal}',
                     f'himidres_{this_gal}']
this_uvh.set_line_products(only=all_line_products)

this_uvh.set_no_cont_products(True)

print(all_line_products)
import time
time.sleep(10)

##############################################################################
# Run staging
##############################################################################


this_uvh.loop_stage_uvdata(do_copy=True, do_contsub=True,
                            do_extract_line=False, do_extract_cont=False,
                            do_remove_staging=False, overwrite=True, strict_config=False)

this_uvh.loop_stage_uvdata(do_copy=False, do_contsub=False,
                            do_extract_line=True, do_extract_cont=False,
                            do_remove_staging=False, overwrite=True, strict_config=False)

this_uvh.loop_stage_uvdata(do_copy=False, do_contsub=False,
                            do_extract_line=False, do_extract_cont=True,
                            do_remove_staging=False, overwrite=True, strict_config=False)

# this_uvh.loop_stage_uvdata(do_copy=False, do_contsub=False,
#                             do_extract_line=False, do_extract_cont=False,
#                             do_remove_staging=True, overwrite=True, strict_config=False)
