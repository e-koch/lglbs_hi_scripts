
'''
Make a series of single channel test images and export a report on the PSF properties.
'''

##############################################################################
# Load routines, initialize handlers
##############################################################################

import os
import sys

from pathlib import Path
import numpy as np


this_galaxy = sys.argv[-3]
this_config = sys.argv[-2]
this_line_product = sys.argv[-1]

# Locate the master key
key_file = "/home/ekoch/lglbs_hi_scripts/lglbs_keys/master_key_fir.txt"


# Set the logging
from phangsPipeline import phangsLogger as pl
pl.setup_logger(level='DEBUG', logfile=None)

# Imports

from phangsPipeline import handlerKeys as kh

from phangsPipeline import handlerTestImaging as tih

# Initialize key handlers
this_kh = kh.KeyHandler(master_key = key_file)
# Make missing directories
this_kh.make_missing_directories(imaging=True, derived=True, postprocess=True, release=True)

##############################################################################
# Set up what we do this run
##############################################################################

this_tih = tih.TestImagingHandler(key_handler=this_kh)

this_tih.set_interf_configs(only=[this_config])
this_tih.set_line_products(only=[this_line_product])


# Define test parameters
robusts = np.round(np.arange(-2.0, 2.1, 0.25), 2)

this_tih.set_test_params(
    weightings=['briggs'] * len(robusts),
    robustnums=robusts,
    tapers_arcsec=None,
)

# Run the test imaging loop
# To remake any images, set overwrite=True
results = this_tih.loop_test_imaging(overwrite=True)

# Generate all diagnostic plots and export CSV in one call
# Exports an HTML report showing the results table and diagnostic plots
# together.
this_tih.make_report()

