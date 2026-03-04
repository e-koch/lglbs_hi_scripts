#!/usr/bin/env python
# 
# Run this script inside CASA!
# 

##############################################################################
# Load routines, initialize handlers
##############################################################################

import os, sys
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
this_imh = imh.ImagingHandler(key_handler = this_kh)
this_pph = pph.PostProcessHandler(key_handler= this_kh)

# Make missing directories
this_kh.make_missing_directories(imaging=True,derived=True,postprocess=True,release=True)

##############################################################################
# Set up what we do this run
##############################################################################

this_targ = 'ic10ctr'

#this_uvh.set_targets(only=['ic10ctr','ic1613ctr','ngc6822','wlmctr'])
this_uvh.set_targets(only=[this_targ])
this_uvh.set_interf_configs(only=['B+C+D','B','C','D'])
this_uvh.set_line_products(only=['oh1612'])
this_uvh.set_no_cont_products(True)

this_imh.set_targets(only=[this_targ])
this_imh.set_interf_configs()
this_imh.set_no_cont_products(True)
this_imh.set_line_products(only=['oh1612'])

this_pph.set_targets()
this_pph.set_interf_configs()
this_pph.set_feather_configs()

# Switches for what steps to run

do_staging = False
do_imaging = True
do_postprocess = False
do_stats = False

##############################################################################
# Run staging
##############################################################################

if do_staging:

    this_uvh.loop_stage_uvdata(do_copy=True, do_contsub=True, 
                               do_extract_line=False, do_extract_cont=False,
                               do_remove_staging=False, overwrite=True, strict_config=False)
    
    this_uvh.loop_stage_uvdata(do_copy=False, do_contsub=False, 
                               do_extract_line=True, do_extract_cont=False,
                               do_remove_staging=False, overwrite=True, strict_config=False)
    
    this_uvh.loop_stage_uvdata(do_copy=False, do_contsub=False, 
                               do_extract_line=False, do_extract_cont=True,
                               do_remove_staging=False, overwrite=True, strict_config=False)
    
    this_uvh.loop_stage_uvdata(do_copy=False, do_contsub=False, 
                               do_extract_line=False, do_extract_cont=False,
                               do_remove_staging=True, overwrite=True, strict_config=False)
    
##############################################################################
# Step through imaging
##############################################################################

# My notes - drop 1000" scale, extend PB cutoff down to 0.5 for 6822
# at least, chanchunks = -1 for CASA 5.7 imaging

if do_imaging:

    this_imh.loop_imaging(do_all = False, do_dirty_image=True)
    #this_imh.loop_imaging(do_all = False, do_dirty_image=True)
    #this_imh.loop_imaging(do_all = True)

##############################################################################
# Step through postprocessing
##############################################################################

if do_postprocess:

    this_pph.loop_postprocess(do_prep = True, do_feather = True, 
                              do_mosaic = True, do_cleanup = True)

