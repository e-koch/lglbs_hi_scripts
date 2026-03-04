'''
Run this script INSIDE CASA or with CASA available.

'''

import os
import sys

####
# Settings
####


this_target = sys.argv[-4]
this_config = sys.argv[-3]
this_line_product = sys.argv[-2]
this_chunksize = int(sys.argv[-1])

print(f"Gathering {this_target} {this_config} {this_line_product} with chunksize {this_chunksize}")

####
# Imports
####

key_file = "/home/erickoch/lglbs_hi_scripts/lglbs_keys/master_key_hydra.txt"

# Set the logging
from phangsPipeline import phangsLogger as pl
pl.setup_logger(level='DEBUG', logfile=None)

from phangsPipeline import handlerKeys as kh
from phangsPipeline.handlerImagingChunked import ImagingChunkedHandler
from phangsPipeline import handlerPostprocess as pph


# Initialize the various handler objects. First initialize the
# KeyHandler, which reads the master key and the files linked in the
# master key. Then feed this keyHandler, which has all the project
# data, into the other handlers (VisHandler, ImagingHandler,
# PostProcessHandler), which run the actual pipeline using the project
# definitions from the KeyHandler.

this_kh = kh.KeyHandler(master_key=key_file)
# Make any missing directories
this_kh.make_missing_directories(imaging=True, derived=True, postprocess=True, release=True)


# Initialize the ImagingChunkedHandler
# Unlike ImageHandler, you must specify the target, configuration and line_name
# ImagingChunkedHandler is designed to image individual products that require
# chunking to process efficiently.
this_imh = ImagingChunkedHandler(this_target,
                                 this_config,
                                 this_line_product,
                                 this_kh,
                                 chunksize=this_chunksize,
                                 make_temp_dir=False)

print(f"Gathering {this_imh.nchunks} chunks")
# When running per chunk, combining into final cubes is a separate call
this_imh.task_complete_gather_into_cubes(root_name='all')

##############################################################################
# Step through postprocessing
##############################################################################

# Postprocess the data in CASA after imaging. This involves primary
# beam correction, linear mosaicking, feathering, conversion to Kelvin
# units, and some downsampling to save space.

this_pph = pph.PostProcessHandler(this_kh)
this_pph.set_targets(only=[this_target])
this_pph.set_interf_configs(only=[this_config])
this_pph.set_line_products(only=[this_line_product])

this_pph.loop_postprocess(do_prep=True, do_feather=True,
                          do_mosaic=True, do_cleanup=True)
