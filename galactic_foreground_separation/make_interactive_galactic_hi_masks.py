
'''
Draw masks around the obvious foreground HI structures in the VLA
cubes.
'''

from pathlib import Path

import numpy as np

import astropy.units as u
from spectral_cube import SpectralCube

from cube_analysis.interactive_masking import make_interactive_cube_mask


# vla_data_path = Path("/reduction10/erickoch/LGLBS/line_imaging/postprocess/")
vla_data_path = Path("/reduction/erickoch/LGLBS/line_imaging/postprocess/")


# galaxy_names = ['ngc6822', 'ic10', 'wlm', 'ic1613']
# hi_cube_name_keys = ['hilores', 'hi21cm_1p2kms', 'himidres']

galaxy_names = ['m31']
hi_cube_name_keys = ['hilores']

# galaxy_names = ['m33']
# hi_cube_name_keys = ['himidres']


for this_gal in galaxy_names:
    print(f"Working on {this_gal}")

    for this_key in hi_cube_name_keys:

        print(f"Working on {this_key}")

        vla_filenames = list((vla_data_path / this_gal).glob(f"{this_gal}*{this_key}_pbcorr_trimmed_k.fits"))

        print(f"Found these files: {[this_name.name for this_name in vla_filenames]}")

        for this_vla_filename in vla_filenames:

            output_maskname = (vla_data_path / this_gal / f"{this_gal}_{this_key}_galacticHI_mask.fits")
            if output_maskname.exists():
                print(f"Skipping {output_maskname}")
                continue

            # Print the velocity range
            print(f"Working on {this_vla_filename}")
            vel_axis = SpectralCube.read(this_vla_filename).with_spectral_unit(u.km/u.s)
            print(f"Velocity range: {np.min(vel_axis.spectral_axis)} to {np.max(vel_axis.spectral_axis)}")

            make_mask = True if input(f"Make mask for {this_vla_filename}? (y/n) ") == 'y' else False

            if not make_mask:
                print(f"Skipping {this_vla_filename}")
                continue

            make_interactive_cube_mask(this_vla_filename, output_maskname,
                                       in_memory=True)

    # CORRECTION: For M31 when smoothed to ~120 pc beam, the mask needs to be expanded a bit to fully capture all foreground HI.

    if this_gal == 'm31':

        # Take channel 29 and expand to 32 with minor erosion to account for M31's rotation.
        with fits.open(output_maskname, mode='update') as hdu:
            hdu[0].data[30:33] = hdu[0].data[29]


