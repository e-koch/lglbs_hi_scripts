# This script demonstrates generation of a "PHANGS product suite" for
# a data set that doesn't have a key file suite set up.

import sys


from spectral_cube import SpectralCube
from radio_beam import Beam
from astropy.io import fits
from astropy.wcs import wcs
import scipy.ndimage as nd
import numpy as np
from astropy import units as u

# Lazy appending to path rather than full install "for now"
# Append to path:
sys.path.append("/home/erickoch/phangs_imaging_scripts/")

from phangsPipeline import scConvolution as scc
from phangsPipeline import scNoiseRoutines as scn
from phangsPipeline import scMaskingRoutines as scm
from phangsPipeline import scMoments as sco

from pathlib import Path

# ####################################################
# Definitions/tuning
# ####################################################


working_dir = Path(sys.argv[-1])

# infile = sys.argv[-2]

outfile = sys.argv[-2]

# dist = "0.84Mpc"

this_base_infile = str(working_dir / outfile)

# Get the beam size
this_beamsize = Beam.from_fits_header(this_base_infile)

# beam_ratios_smooth = [1.5, 2., 4.]
# target_beams = np.array(beam_ratios_smooth) * this_beamsize.major.to(u.arcsec).value

# target_angular = [f'{this_targbeam:.1f}arcsec' for this_targbeam in target_beams]
target_angular = []
target_physical = []

print(f"Running on angular scales: {target_angular}")

do_noise = False
do_strictmask = False
do_broadmask = True
do_moments = True

# target_angular = ['8arcsec','15arcsec','21arcsec']
# target_physical = ['150pc','500pc','1000pc','1500pc']

mask_kwargs = {'hi_thresh': 4.0, 'lo_thresh': 2.0,
               'hi_nchan': 2, 'lo_nchan': 2,
               'min_beams': 0.75}


mom_list = {
    'strictmom0':{
        'algorithm':'mom0',
        'mask':'_strictmask.fits',
        'ext':'_strict_mom0.fits',
        'ext_error':'_strict_emom0.fits',
        'kwargs':{},
        },
    'broadmom0':{
        'algorithm':'mom0',
        'mask':'_broadmask.fits',
        'ext':'_broad_mom0.fits',
        'ext_error':'_broad_emom0.fits',
        'kwargs':{},
        },
    'strictmom1':{
        'algorithm':'mom1',
        'mask':'_strictmask.fits',
        'ext':'_strict_mom1.fits',
        'ext_error':'_strict_emom1.fits',
        'kwargs':{},
        },
    'strictmom2':{
        'algorithm':'mom2',
        'mask':'_strictmask.fits',
        'ext':'_strict_mom2.fits',
        'ext_error':'_strict_emom2.fits',
        'kwargs':{},
        },
    'strictew':{
        'algorithm':'ew',
        'mask':'_strictmask.fits',
        'ext':'_strict_ew.fits',
        'ext_error':'_strict_eew.fits',
        'kwargs':{},
        },
    'tpeak':{
        'algorithm':'tpeak',
        'mask':'_broadmask.fits',
        'ext':'_broad_tpeak.fits',
        'ext_error':'_broad_etpeak.fits',
        'kwargs':{},
        },
    }


this_infile = working_dir / outfile
this_base_infile = str(working_dir / outfile)

this_base_cube = SpectralCube.read(this_base_infile)
this_base_cube.allow_huge_operations = True


# ####################################################
# Noise Estimation
# ####################################################

if do_noise:

    # native
    print("Noise native res for ", this_base_infile)
    scn.recipe_phangs_noise(
        incube=this_base_cube,
        outfile=this_base_infile.replace('.fits','_noise.fits'),
        overwrite=True)

# ####################################################
# Strict Mask Creation
# ####################################################

# this_base_infile = str(working_dir / outfile)

if do_strictmask:

    import time

    t0 = time.time()

    # native
    this_infile = this_base_infile
    print("Strict mask for ", this_infile)
    scm.recipe_phangs_strict_mask(
        this_base_cube,
        this_infile.replace('.fits','_noise.fits'),
        outfile=this_infile.replace('.fits','_strictmask.fits'),
        mask_kwargs=mask_kwargs,
        overwrite=True)

    # Now filter out Galactic HI emission that remains
    # in the VLA cube (particularly for NGC6822)

    this_mask_file = this_infile.replace('.fits','_strictmask.fits')

    t1 = time.time()

    print(f"Strict masking took: {t1-t0}")

# ####################################################
# Broad Mask Creation
# ####################################################

if do_broadmask:

    this_base_infile = str(working_dir / outfile)
    template_mask = this_base_infile.replace('.fits','_strictmask.fits')
    broad_mask_file = this_base_infile.replace('.fits','_broadmask.fits')

    list_of_masks = [template_mask]

    for this_ext in target_angular:
        this_infile = this_base_infile.replace('.fits', f'_{this_ext}.fits')
        this_mask = this_infile.replace('.fits','_strictmask.fits')
        list_of_masks.append(this_mask)

    for this_ext in target_physical:
        this_infile = this_base_infile.replace('.fits', f'_{this_ext}.fits')
        this_mask = this_infile.replace('.fits','_strictmask.fits')
        list_of_masks.append(this_mask)

    scm.recipe_phangs_broad_mask(
        template_mask,
        outfile=broad_mask_file,
        list_of_masks = list_of_masks,
        grow_xy = 5,
        grow_v = 4,
        return_spectral_cube=False, overwrite=True,
        recipe='anyscale', fraction_of_scales=0.25)

# ####################################################
# Moments
# ####################################################

if do_moments:

    # native
    this_base_infile = str(working_dir / outfile)
    this_infile = this_base_infile
    print("Moments for ", this_infile)

    for this_mom_key in mom_list:
        this_mom = mom_list[this_mom_key]
        if this_mom['mask'] == '_broadmask.fits':
            this_mask_file = this_base_infile.replace('.fits',this_mom['mask'])
        else:
            this_mask_file = this_infile.replace('.fits',this_mom['mask'])
        this_noise_file = this_infile.replace('.fits','_noise.fits')
        this_out_file = this_infile.replace('.fits',this_mom['ext'])
        this_error_file = this_infile.replace('.fits',this_mom['ext_error'])

        sco.moment_generator(
            this_infile,
            mask=this_mask_file, noise=this_noise_file,
            moment=this_mom['algorithm'], momkwargs=this_mom['kwargs'],
            outfile=this_out_file, errorfile=this_error_file,
            channel_correlation=None,
        )
