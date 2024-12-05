
import sys
import os

from astropy.io import fits

this_target = sys.argv[-3]
this_config = sys.argv[-2]
this_line_product = sys.argv[-1]

from phangsPipeline import phangsLogger as pl
pl.setup_logger(level='DEBUG', logfile=None)


from phangsPipeline import handlerKeys as kh
from phangsPipeline import handlerDerived as der

# Initialize key handler
key_file = "/home/erickoch/lglbs_hi_scripts/lglbs_keys/master_key_hydra.txt"

this_kh = kh.KeyHandler(master_key=key_file)

this_der = der.DerivedHandler(key_handler=this_kh)


this_der.set_targets(only=[this_target])
# this_der.set_interf_configs(only=['C','C+D'])
this_der.set_interf_configs()
this_der.set_line_products(only=[this_line_product])
# this_der.set_line_products(only=['hilores'])
# this_der.set_feather_configs(only=[f"{this_config}+tp"])
this_der.set_no_cont_products(True)

# Set steps
do_convolve = True
do_noise = True
do_strictmask = True
do_broadmask = True
do_moments = True
do_secondary = True

# this_config = f"{this_config}+tp"

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
        mask_filename = f"{this_kh._postprocess_root}/{this_target}/{this_target}_{this_line}_galacticHI_mask.fits"
        print(mask_filename)

        if not os.path.exists(mask_filename):
            # raise Exception(f"Missing {mask_filename}")
            print(f"Missing {mask_filename}. Skipping.")
            continue

        # Otherwise open the mask and exclude any foreground components from the strictmask
        mask_fg = fits.open(mask_filename)[0].data > 0

        # Native resolution
        strictmask_filename = f"{this_kh._derived_root}/{this_target}/{this_target}_{this_config}_{this_line}_strictmask.fits"
        if not os.path.exists(strictmask_filename):
            raise Exception(f"Missing {strictmask_filename}")
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

if do_moments:
    this_der.loop_derive_products(do_convolve = False, do_noise = False,
                                  do_strictmask = False, do_broadmask = False,
                                  do_moments = True, do_secondary = False)

if do_secondary:
    this_der.loop_derive_products(do_convolve = False, do_noise = False,
                                  do_strictmask = False, do_broadmask = False,
                                  do_moments = False, do_secondary = True)
