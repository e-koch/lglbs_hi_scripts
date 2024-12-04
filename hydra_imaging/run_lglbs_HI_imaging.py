
'''
Make dirty images of the wide-field C+D config data for the dwarfs.

Produces dirty images for:
- 4 OH lines
- HI at 4.2 km/s (10x native)
- HI at 2.1 km/s (5x native)
- HI at 0.84 km/s (2x native)
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

from pathlib import Path
import numpy as np



this_idx = int(sys.argv[-1])
this_chunksize = int(sys.argv[-2])
this_config = sys.argv[-3]
this_line_product = sys.argv[-4]
this_galaxy = sys.argv[-5]

# -1 offset from job array index
this_chunk_num = this_idx-1

# Locate the master key
key_file = "/home/erickoch/lglbs_hi_scripts/lglbs_keys/master_key_hydra.txt"


# Set the logging
from phangsPipeline import phangsLogger as pl
pl.setup_logger(level='DEBUG', logfile=None)

# Imports

#sys.path.insert(1, )
from phangsPipeline import handlerKeys as kh

from phangsPipeline.handlerImagingChunked import ImagingChunkedHandler

# Initialize key handlers
this_kh = kh.KeyHandler(master_key = key_file)
# Make missing directories
this_kh.make_missing_directories(imaging=True, derived=True, postprocess=True, release=True)

##############################################################################
# Set up what we do this run
##############################################################################

this_imh = ImagingChunkedHandler(this_galaxy,
                                 this_config,
                                 this_line_product,
                                 this_kh,
                                 chunksize=this_chunksize,
                                 make_temp_dir=True)


print(f"Chunk {this_chunk_num} of {this_imh.nchunks}")

this_imh.run_imaging(chunk_num=this_chunk_num,
                     do_all=True,)


