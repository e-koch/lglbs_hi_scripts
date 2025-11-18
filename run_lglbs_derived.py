
##############################################################################
# Load routines, initialize handlers
##############################################################################

import os
import sys

from astropy.io import fits
import numpy as np

# Pass the target name from the cmd line

this_target = sys.argv[-1]

this_config = 'C+D+tp'

all_targets = ['ngc6822','wlm','ic10','ic1613', 'm33', 'm31']
# all_targets = ['m31']

if this_target not in all_targets:
    raise Exception(f"Please specify a valid target. Must be one of {all_targets}.")

# Script directory
script_dir = "/home/erickoch/lglbs_hi_scripts/"

# Locate the master key
# key_file = '/data/tycho/0/leroy.42/reduction/vla/lglbs_pipeline_configs/lglbs_keys/master_key.txt'
key_file = f'{script_dir}/lglbs_keys/master_key.txt'

sys.path.append(os.path.expanduser("~/phangs_imaging_scripts"))

# Set the logging
from phangsPipeline import phangsLogger as pl
pl.setup_logger(level='DEBUG', logfile=None)

# Imports

#sys.path.insert(1, )
from phangsPipeline import handlerKeys as kh
from phangsPipeline import handlerDerived as der

# Initialize key handler

this_kh = kh.KeyHandler(master_key=key_file)
this_der = der.DerivedHandler(key_handler=this_kh)

# Make missing directories

this_kh.make_missing_directories(imaging=True,
                                 derived=True,
                                 postprocess=True,
                                 release=True)

##############################################################################
# Set up what we do this run
##############################################################################

#this_der.set_targets(only=['ngc6822','wlm','ic10'])
this_der.set_targets(only=[this_target])
# this_der.set_interf_configs(only=['C','C+D'])
this_der.set_interf_configs()
# this_der.set_line_products(only=['hilores', 'hi21cm_1p2kms', 'himidres'])
# this_der.set_line_products(only=['hi21cm_0p8kms'])
this_der.set_line_products(only=['hilores'])
# this_der.set_line_products(only=['himidres'])
this_der.set_feather_configs(only=[this_config])
this_der.set_no_cont_products(True)

# Set steps
do_convolve = False
do_noise = False
do_strictmask = True
do_broadmask = True
do_moments = True
do_secondary = True


##############################################################################
# Step through derived product creation
##############################################################################

if do_convolve:
    this_der.loop_derive_products(do_convolve = True, do_noise = False,
                                  do_strictmask = False, do_broadmask = False,
                                  do_moments = False, do_secondary = False)

if do_noise:
    this_der.loop_derive_products(do_convolve = False, do_noise = True,
                                  do_strictmask = False, do_broadmask = False,
                                  do_moments = False, do_secondary = False)

if do_strictmask:
    this_der.loop_derive_products(do_convolve = False, do_noise = False,
                                  do_strictmask = True, do_broadmask = False,
                                  do_moments = False, do_secondary = False)

    # Check for the existence of a MW HI foreground mask
    for this_line in this_der.get_line_products():

        # Native resolution
        strictmask_filename = f"{this_kh._derived_root}/{this_target}/{this_target}_{this_config}_{this_line}_strictmask.fits"
        if not os.path.exists(strictmask_filename):
            raise Exception(f"Missing {strictmask_filename}")


        # Read in any custom regions to exclude
        reg_filename = f"{script_dir}/{this_target}_artifact_mask.reg"
        if os.path.exists(reg_filename):
            print(f"Using {reg_filename}")

            from regions import Regions
            from astropy.wcs import WCS
            regs = Regions.read(reg_filename, format='ds9')

            mask_reg = np.zeros(fits.open(strictmask_filename)[0].shape[1:], dtype=bool)
            for this_reg in regs:
                reg_pix = this_reg.to_pixel(WCS(fits.getheader(strictmask_filename)).celestial)
                mask_reg = np.logical_or(mask_reg,
                                         reg_pix.to_mask().to_image(fits.open(strictmask_filename)[0].shape[1:]))
        else:
            regs = None


        mask_filename = f"{this_kh._postprocess_root}/{this_target}/{this_target}_{this_line}_galacticHI_mask.fits"
        print(mask_filename)

        if not os.path.exists(mask_filename) and regs is None:
            # raise Exception(f"Missing {mask_filename}")
            print(f"Missing {mask_filename}. Skipping.")
            continue

        # Otherwise open the mask and exclude any foreground components from the strictmask
        mask_fg = fits.open(mask_filename)[0].data > 0

        if regs is not None:
            for ii in range(mask_fg.shape[0]):
                mask_fg[ii][mask_reg] = True

        # Open the mask
        with fits.open(strictmask_filename, 'update') as hdul:
            hdul[0].data[mask_fg] = False
            hdul.flush()

        res_dict = this_kh.get_ang_res_dict(config=this_config, product=this_line)
        res_list = list(res_dict)
        if len(res_list) > 0:
            res_list.sort()
        for this_res_tag in res_list:
            strictmask_filename = f"{this_kh._derived_root}/{this_target}/{this_target}_{this_config}_{this_line}_{this_res_tag}_strictmask.fits"

            if not os.path.exists(strictmask_filename):
                raise Exception(f"Missing {strictmask_filename}")
            # Open the mask
            with fits.open(strictmask_filename, 'update') as hdul:
                hdul[0].data[mask_fg] = False
                hdul.flush()

        res_dict = this_kh.get_phys_res_dict(config=this_config, product=this_line)
        res_list = list(res_dict)
        if len(res_list) > 0:
            res_list.sort()
        for this_res_tag in res_list:
            strictmask_filename = f"{this_kh._derived_root}/{this_target}/{this_target}_{this_config}_{this_line}_{this_res_tag}_strictmask.fits"

            if not os.path.exists(strictmask_filename):
                raise Exception(f"Missing {strictmask_filename}")
            # Open the mask
            with fits.open(strictmask_filename, 'update') as hdul:
                hdul[0].data[mask_fg] = False
                hdul.flush()

if do_broadmask:
    this_der.loop_derive_products(do_convolve = False, do_noise = False,
                                  do_strictmask = False, do_broadmask = True,
                                  do_moments = False, do_secondary = False)

    # Enforce MW foreground mask in all cases:

    # Check for the existence of a MW HI foreground mask
    for this_line in this_der.get_line_products():

        mask_filename = f"{this_kh._postprocess_root}/{this_target}/{this_target}_{this_line}_galacticHI_mask.fits"
        print(mask_filename)

        if not os.path.exists(mask_filename) and regs is None:
            # raise Exception(f"Missing {mask_filename}")
            print(f"Missing {mask_filename}. Skipping.")
            continue

        # Otherwise open the mask and exclude any foreground components from the strictmask
        mask_fg = fits.open(mask_filename)[0].data > 0

        # Native resolution
        broadmask_filenames = []
        broadmask_filename = f"{this_kh._derived_root}/{this_target}/{this_target}_{this_config}_{this_line}_broadmask.fits"
        if not os.path.exists(broadmask_filename):
            raise Exception(f"Missing {broadmask_filename}")
        broadmask_filenames.append(broadmask_filename)

        res_dict = this_kh.get_phys_res_dict(config=this_config, product=this_line)
        res_list = list(res_dict)
        if len(res_list) > 0:
            res_list.sort()
        res_dict_ang = this_kh.get_ang_res_dict(config=this_config, product=this_line)
        res_list_ang = list(res_dict_ang)
        if len(res_list_ang) > 0:
            res_list_ang.sort()
        res_list.extend(res_list_ang)

        for this_res_tag in res_list:
            broadmask_filename = f"{this_kh._derived_root}/{this_target}/{this_target}_{this_config}_{this_line}_{this_res_tag}_broadmask.fits"

            if os.path.exists(broadmask_filename):
                broadmask_filenames.append(broadmask_filename)

        # Loop through all filenames
        for broadmask_filename in broadmask_filenames:
            # Open the mask
            with fits.open(broadmask_filename, 'update') as hdul:
                hdul[0].data[mask_fg] = False
                hdul.flush()



if do_moments:
    this_der.loop_derive_products(do_convolve = False, do_noise = False,
                                  do_strictmask = False, do_broadmask = False,
                                  do_moments = True, do_secondary = False)

if do_secondary:
    this_der.loop_derive_products(do_convolve = False, do_noise = False,
                                  do_strictmask = False, do_broadmask = False,
                                  do_moments = False, do_secondary = True)
