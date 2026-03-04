
'''
Adds manual masking for small continuum subtraction residuals towards bright continuum sources
And subtracts the "residual" continuum pedestal within the frequency range the rflag was disabled

'''

import warnings
import os
from pathlib import Path

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import stats

from astropy.io import fits
from astropy.wcs.utils import proj_plane_pixel_scales

from astropy.modeling.models import Gaussian2D
from astropy.modeling import models, fitting


from regions import EllipseSkyRegion

from spectral_cube import SpectralCube

import numpy as np
import matplotlib.pyplot as plt


def fit_linear_background(subcube, plot=False):

    spectral_axis = subcube.spectral_axis.value

    # Sigma clip to remove signal, then take the median
    data_clipped = stats.sigma_clip(subcube.unitless_filled_data[:],
                                    axis=0,
                                    sigma=2.0,
                                    maxiters=5,
                                    cenfunc='median',
                                    stdfunc='mad_std')

    med_clipped = np.nanmedian(data_clipped, axis=0)

    out_model = np.zeros(subcube.shape)
    out_params = np.zeros((2,) + subcube.shape[1:])

    for y, x in np.ndindex(subcube.shape[1], subcube.shape[2]):
        if np.isnan(data_clipped[:, y, x]).all():
            continue

        line_orig = models.Linear1D(slope=0.0,
                                    intercept=med_clipped[y, x])

        fitter = fitting.LinearLSQFitter()

        mask = np.isfinite(data_clipped[:, y, x])
        line_fit = fitter(line_orig,
                          spectral_axis[mask],
                          data_clipped[:, y, x][mask])

        out_model[:, y, x] = line_fit(spectral_axis)
        out_params[:, y, x] = line_fit.parameters

    # if plot:
    #     fig, axes = plt.subplots(subcube.shape[1],
    #                              subcube.shape[2],
    #                              figsize=(15, 15),
    #                              sharex=True,
    #                              sharey=True)

    #     for y, x in np.ndindex(subcube.shape[1], subcube.shape[2]):
    #         axes[y, x].plot(spectral_axis, subcube[:, y, x].value)
    #         axes[y, x].plot(spectral_axis, out_model[:, y, x])


    return out_model, out_params


def fit_gaussian_background(subcube):

    # Sigma clip to remove signal, then take the median
    data_clipped = stats.sigma_clip(subcube.unitless_filled_data[:],
                                    axis=0,
                                    sigma=2.0,
                                    maxiters=5,
                                    cenfunc='median',
                                    stdfunc='mad_std')

    med_clipped = np.nanmedian(data_clipped, axis=0)
    med_clipped = np.nan_to_num(med_clipped, copy=False)

    # Define the spatial grid for the fit centered at y, x = 32, 32
    yy, xx = np.mgrid[:subcube.shape[1], :subcube.shape[2]]
    size = subcube.shape[1] // 2

    # Define a single 2D Gaussian model.
    p_init_gauss2D = models.Gaussian2D(x_mean=xx[size, size],
                                       y_mean=yy[size, size],
                                       amplitude=med_clipped.max(),
                                       x_stddev=6.0,
                                       y_stddev=6.0)

    # And fit with the Levenberg-Marquardt algorithm and least squares statistic.
    fit_p = fitting.LevMarLSQFitter()

    with warnings.catch_warnings():
        # Ignore model linearity warning from the fitter
        warnings.simplefilter('ignore')
        p_gauss2D = fit_p(p_init_gauss2D, xx, yy, med_clipped)

    model = p_gauss2D(xx, yy)

    return model



data_path = Path("/reduction/erickoch/LGLBS/line_imaging/postprocess/")

# All in km/s
vel_ranges = {"m31": [60, -650] * u.km / u.s,   # [50, -625]
              "m33": [40, -340] * u.km / u.s,
              "ngc6822": [50, -625] * u.km / u.s,
              "wlm": [-50, -210] * u.km / u.s,
              "ic10": [-240, -450] * u.km / u.s,
              "ic1613": [-170, -290] * u.km / u.s,}


# No strong residual continuum sources in:
# IC10
# IC1613
# NGC6822


# M33

this_gal = "m33"
this_gal_datapath = data_path / this_gal

cont_coords = SkyCoord(ra=["1h32m22.8"],
                       dec=["+30d44m06.2"],
                       frame='fk5',
                       unit=(u.hourangle, u.deg))

# Loop over all "_pbcorr_trimmed_k.fits" files in the postprocess directory
# for...

this_filename = this_gal_datapath / "m33_C+D+tp_himidres_pbcorr_trimmed_k.fits"

# Make a backup before we fiddle with the data.
this_filename_backup = Path(f"{this_filename}_backup")

if not this_filename_backup.exists():
    print(f"Backing up {this_filename}")
    os.system(f"cp {this_filename} {this_filename_backup}")
else:
    print(f"Already backed up {this_filename}")


cube = SpectralCube.read(this_filename)

nbeams = 5

# For loop through all continuum coordinates
for this_coord in cont_coords:
    this_region = EllipseSkyRegion(center=this_coord,
                                   width=nbeams * cube.beam.major.to(u.deg),
                                   height=nbeams * cube.beam.major.to(u.deg),
                                   angle=cube.beam.pa)

    this_pix_region = this_region.to_pixel(cube.wcs.celestial)
    this_mask = this_pix_region.to_mask()

    slices_large, slices_small = this_mask.get_overlap_slices(cube.shape[1:])

    subcube = cube[(slice(None),) + slices_large].with_mask(this_mask.data.astype(bool))

    # Cut to the velocity range
    # Pad by ~2 channels on each side
    low_chan = max(subcube.closest_spectral_channel(vel_ranges[this_gal][0]) - 1, 0)
    high_chan = min(subcube.closest_spectral_channel(vel_ranges[this_gal][1]) + 1, cube.shape[0]-1)
    subcube = subcube[low_chan:high_chan]


    # model = fit_gaussian_background(subcube)

    # Plot the fit and residual
    # fig, axs = plt.subplots(1, 3, figsize=(8, 4),
    #                         sharex=True, sharey=True)
    # im = axs[0].imshow(med_clipped, origin='lower')
    # fig.colorbar(im, ax=axs[0])
    # im2 = axs[1].imshow(model, origin='lower')
    # fig.colorbar(im2, ax=axs[1])
    # im3 = axs[2].imshow(med_clipped - model, origin='lower')
    # fig.colorbar(im3, ax=axs[2])


    model, params = fit_linear_background(subcube)

    # Open the FITS cube and make + update the data.
    with fits.open(this_filename, 'update') as hdul:
        hdul[0].data[(slice(low_chan, high_chan),) + slices_large] -= model


# M31

cont_coords = SkyCoord(ra=["0h42m18.9",
                           "0h43m55.0",
                           "0h42m17.0",
                           "0h37m37.0",
                           "0h38m48.4",
                           "0h46m48.3",],
                       dec=["+41d29m26.9",
                            "+40d46m34.0",
                            "+40d09m59.9",
                            "+39d38m11.0",
                            "+41d16m07.8",
                            "+42d08m56.9"],
                       frame='fk5',
                       unit=(u.hourangle, u.deg))

this_gal = "m31"

this_gal_datapath = data_path / this_gal

this_filename = this_gal_datapath / "m31_C+D+tp_hilores_pbcorr_trimmed_k.fits"

# Make a backup before we fiddle with the data.
this_filename_backup = Path(f"{this_filename}_backup")

if not this_filename_backup.exists():
    print(f"Backing up {this_filename}")
    os.system(f"cp {this_filename} {this_filename_backup}")
else:
    print(f"Already backed up {this_filename}")


cube = SpectralCube.read(this_filename)

nbeams = 5

# For loop through all continuum coordinates
for this_coord in cont_coords:
    this_region = EllipseSkyRegion(center=this_coord,
                                   width=nbeams * cube.beam.major.to(u.deg),
                                   height=nbeams * cube.beam.major.to(u.deg),
                                   angle=cube.beam.pa)

    this_pix_region = this_region.to_pixel(cube.wcs.celestial)
    this_mask = this_pix_region.to_mask()

    slices_large, slices_small = this_mask.get_overlap_slices(cube.shape[1:])

    subcube = cube[(slice(None),) + slices_large].with_mask(this_mask.data.astype(bool))

    # Cut to the velocity range
    # Pad by ~2 channels on each side
    low_chan = max(subcube.closest_spectral_channel(vel_ranges[this_gal][0]) - 1, 0)
    high_chan = min(subcube.closest_spectral_channel(vel_ranges[this_gal][1]) + 1, cube.shape[0]-1)
    subcube = subcube[low_chan:high_chan]


    model, params = fit_linear_background(subcube, plot=False)

    # input("Press enter to continue")
    # plt.close()

    # Open the FITS cube and make + update the data.
    with fits.open(this_filename, 'update') as hdul:
        hdul[0].data[(slice(low_chan, high_chan),) + slices_large] -= model


# WLM


this_gal = "wlm"
this_gal_datapath = data_path / this_gal

cont_coords = SkyCoord(ra=["0h01m41.9",
                           "0h3m27.3",
                           "0h03m16.1",],
                       dec=["-15d40m38.1",
                            "-15d47m06.4",
                            "-15d30m13.2"],
                       frame='fk5',
                       unit=(u.hourangle, u.deg))

# Loop over all "_pbcorr_trimmed_k.fits" files in the postprocess directory
# for...

for this_filename in this_gal_datapath.glob("*_pbcorr_trimmed_k.fits"):

    # Make a backup before we fiddle with the data.
    this_filename_backup = Path(f"{this_filename}_backup")

    if not this_filename_backup.exists():
        print(f"Backing up {this_filename}")
        os.system(f"cp {this_filename} {this_filename_backup}")
    else:
        print(f"Already backed up {this_filename}")


    cube = SpectralCube.read(this_filename)

    nbeams = 5


    # For loop through all continuum coordinates
    for this_coord in cont_coords:
        this_region = EllipseSkyRegion(center=this_coord,
                                    width=nbeams * cube.beam.major.to(u.deg),
                                    height=nbeams * cube.beam.major.to(u.deg),
                                    angle=cube.beam.pa)

        this_pix_region = this_region.to_pixel(cube.wcs.celestial)
        this_mask = this_pix_region.to_mask()

        slices_large, slices_small = this_mask.get_overlap_slices(cube.shape[1:])

        subcube = cube[(slice(None),) + slices_large].with_mask(this_mask.data.astype(bool))

        # Cut to the velocity range
        # Pad by ~2 channels on each side
        low_chan = max(subcube.closest_spectral_channel(vel_ranges[this_gal][0]) - 1, 0)
        high_chan = min(subcube.closest_spectral_channel(vel_ranges[this_gal][1]) + 1, cube.shape[0]-1)
        subcube = subcube[low_chan:high_chan]


        model, params = fit_linear_background(subcube, plot=False)

        # input("Press enter to continue")
        # plt.close()

        # Open the FITS cube and make + update the data.
        with fits.open(this_filename, 'update') as hdul:
            hdul[0].data[(slice(low_chan, high_chan),) + slices_large] -= model
