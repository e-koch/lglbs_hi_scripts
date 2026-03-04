
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

all_targs = ['ic10ctr', 'ic1613ctr', 'ngc6822', 'wlmctr']

if not this_targ in all_targs:
    raise ValueError(f"Cannot find target {this_targ} in list of target names: {all_targs}")

#this_uvh.set_targets(only=['ic10ctr','ic1613ctr','ngc6822','wlmctr'])
this_uvh.set_targets(only=[this_targ])
this_uvh.set_interf_configs(only=['A', 'B', 'C', 'D', 'B+C+D', 'A+B+C+D', 'A+B+C'])

all_line_products = ['hi21cm_1p2kms', 'hi21cm', 'hilores', 'himidres',
                     'oh1612', 'oh1665', 'oh1667', 'oh1720',]

# this_uvh.set_line_products(only=['hilores'])

this_uvh.set_no_cont_products(True)

##############################################################################
# Run staging
##############################################################################

for this_line in all_line_products:

    this_uvh.set_line_products(only=[this_line])

    this_uvh.loop_stage_uvdata(do_copy=True, do_contsub=True,
                                do_extract_line=False, do_extract_cont=False,
                                do_remove_staging=False, overwrite=True, strict_config=False)

    this_uvh.loop_stage_uvdata(do_copy=False, do_contsub=False,
                                do_extract_line=True, do_extract_cont=False,
                                do_remove_staging=False, overwrite=True, strict_config=False)

    this_uvh.loop_stage_uvdata(do_copy=False, do_contsub=False,
                                do_extract_line=False, do_extract_cont=False,
                                do_remove_staging=True, overwrite=True, strict_config=False)
