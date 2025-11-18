
import os
from pathlib import Path
import warnings

import numpy as np

from astropy.io import fits
import astropy.units as u
from astropy.convolution import Gaussian1DKernel

from spectral_cube import SpectralCube
from radio_beam import Beam

from uvcombine import feather_simple_cube

from tqdm import tqdm

# sd_data_path = Path("/reduction10/erickoch/LGLBS/hi_feathering/")
# vla_data_path = Path("/reduction10/erickoch/LGLBS/C+D_HI_2023/")
sd_data_path = Path("/reduction/erickoch/LGLBS/hi_feathering/")
# vla_data_path = Path("/reduction/erickoch/LGLBS/line_imaging/imaging/m33/")
vla_data_path = Path("/reduction/erickoch/LGLBS/line_imaging/imaging/m31/")


galaxy_dict = {\
               # 'm33': ['m33_gbt_vlsr.fits', 1.10],
               'm31': ["M31_GBT_corrected.fits", 1.0],
               # 'ngc6822': ["NGC6822-center_cube_ircs_no_galactic_component.fits",
               #            0.92],
               # 'ic10': ['IC10_GBT_Jy.fits',
               #          0.96],
               # 'wlm': ['WLM_GBT.FITS',
               #         1.14],
               # 'ic1613': ['IC1613_GBT_vegas_K_noresample_lsrk.fits',
               #            0.89],
               }

# Loop through by keys by VLA name. This allows us to run the HI cubes
# made at different spectral resolution
hi_cube_name_keys = ['hilores', 'hi21cm_1p2kms', 'himidres']
# hi_cube_name_keys = ['hi21cm_1p2kms', 'himidres']

for this_key in hi_cube_name_keys:

    # vla_filenames = list(vla_data_path.glob(f"*{this_key}.fits"))
    vla_filenames = list(vla_data_path.glob(f"*{this_key}.fits"))
    print(f"Found these files: {[this_name.name for this_name in vla_filenames]}")

    # for this_gal in galaxy_dict:
    for this_vla_filename in vla_filenames:


        this_gal = this_vla_filename.name.split("_")[0]

        print(f"On {this_gal}: {this_vla_filename.name}")

        if this_gal not in galaxy_dict:
            print(f"Unable to find matching key for {this_vla_filename.name}. Skipping.")
            continue

        this_gbt_filename = galaxy_dict[this_gal][0]
        this_scfactor = galaxy_dict[this_gal][1]

        # this_vla_filename = galaxy_dict[this_gal][1]
        # this_scfactor = galaxy_dict[this_gal][2]

        #this_specslice_low, this_specslice_high = galaxy_dict[this_gal][2]

        vla_cube = SpectralCube.read(vla_data_path / this_vla_filename,
                                     use_dask=True)
        vla_cube.allow_huge_operations = True

        # VLA spatial mask. Account for coverage across all channels (in case some are empty)
        vla_spatial_mask = np.any(vla_cube.mask.include(), axis=0)

        # Nick's reprocessed and gridded GBT cube.
        gbt_cube = SpectralCube.read(sd_data_path / this_gbt_filename,
                                     use_dask=True)
        gbt_cube.allow_huge_operations = True
        print(f"GBT original bunit: {gbt_cube.unit}")

        # Use the proper beam model size, not the one in the header!
        if this_gal != "m33":

            gbt_beam_model = Beam(area=3.69e5 *u.arcsec**2)
            gbt_beam_model.major.to(u.arcmin)

            print(f"Using proper GBT beam area model with {gbt_beam_model.major.to(u.arcmin)} arcmin FWHM")
        else:
            print(f"Using M33 GBT beam model with 9.8 arcmin FWHM")
            # See Koch+18 for difference in the gridding kernel used for M33
            gbt_beam_model = Beam(major=9.8*u.arcmin)

        vel_unit = u.m / u.s

        gbt_cube = gbt_cube.with_beam(gbt_beam_model, raise_error_jybm=False)
        gbt_cube = gbt_cube.with_spectral_unit(vel_unit, velocity_convention='radio')

        # There is an issue with the F2F optical to radio conversion that is not being
        # handled during reprojection. Just set to VRAD if that's the case.
        if "F2F" in gbt_cube.wcs.wcs.ctype[2]:
            gbt_cube.wcs.wcs.ctype[2] = "VRAD"
            gbt_cube.mask._wcs.wcs.ctype[2] = "VRAD"


        fwhm_factor = np.sqrt(8*np.log(2))
        current_resolution = np.abs(np.diff(gbt_cube.spectral_axis)[0]).to(vel_unit)
        target_resolution = np.abs(np.diff(vla_cube.spectral_axis)[0]).to(vel_unit)

        if current_resolution < target_resolution:
            gaussian_width = ((target_resolution**2 - current_resolution**2)**0.5 /
                            current_resolution / fwhm_factor)
            print(f"Spectral smoothing by: {gaussian_width}")
            kernel = Gaussian1DKernel(gaussian_width.value)
            gbt_cube_specsmooth = gbt_cube.spectral_smooth(kernel)
        else:
            gbt_cube_specsmooth = gbt_cube


        gbt_cube_specinterp = gbt_cube_specsmooth.spectral_interpolate(vla_cube.spectral_axis)

        # This was for sanity checking issues that were happening with the reprojection (see below)
        # gbt_cube_specinterp.to(u.K).write(sd_data_path / f"{this_gbt_filename[:-5]}_vlamatch_{this_key}_specinterp.fits",
        #                                           overwrite=True)

        # Grid the GBT data to the VLA grid
        ## The full cube reprojection is failing in some cases. I don't understand why.
        # target_hdr = vla_cube.header.copy()
        # target_hdr['TIMESYS'] = target_hdr['TIMESYS'].lower()
        # gbt_cube_specinterp_reproj = gbt_cube_specinterp.reproject(target_hdr)
        # gbt_cube_specinterp_reproj.to(u.K).write(sd_data_path / f"{this_gbt_filename[:-5]}_vlamatch_{this_key}_reproj.fits",
        #                                           overwrite=True)

        # Do a per-channel version to avoid the problem
        gbt_reproj_filename = sd_data_path / f"{this_gbt_filename[:-5]}_vlamatch_{this_key}.fits"
        # Generate a copy of the VLA cube
        if gbt_reproj_filename.exists():
            gbt_reproj_filename.unlink()

        print(f"Copying {vla_data_path / this_vla_filename} to {gbt_reproj_filename}")
        os.system(f"cp {vla_data_path / this_vla_filename} {gbt_reproj_filename}")

        # Check the units are correct for the GBT data:
        gbt_bunit = gbt_cube.unit
        print(f"GBT original bunit: {gbt_bunit}")
        with fits.open(gbt_reproj_filename, mode="update") as hdulist:
            hdulist[0].header['BUNIT'] = str(gbt_bunit)
            hdulist[0].header.update(gbt_cube.beam.to_header_keywords())

            if len(hdulist) > 1:
                # Remove the beam table
                del hdulist[1]

            hdulist.flush()

        print("Per channel reprojection")
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="WCS1 is missing card")
            for this_chan in tqdm(range(vla_cube.shape[0])):
                reproj_chan = gbt_cube_specinterp[this_chan].reproject(vla_cube[this_chan].header)

                # Apply the VLA spatial mask:
                reproj_chan[~vla_spatial_mask] = np.nan

                with fits.open(gbt_reproj_filename, mode="update") as hdulist:
                    hdulist[0].data[this_chan] = reproj_chan

                    hdulist.flush()

        # Allow reading in the whole cube.
        gbt_cube_specinterp_reproj = SpectralCube.read(gbt_reproj_filename,
                                                       use_dask=True,)
        gbt_cube_specinterp_reproj.allow_huge_operations = True

        # And specifically apply the same PB coverage
        # gbt_cube_specinterp_reproj = gbt_cube_specinterp_reproj.with_mask(vla_spatial_mask)

        # gbt_cube_specinterp_reproj.to(u.K).write(sd_data_path / f"{this_gbt_filename[:-5]}_vlamatch_{this_key}.fits",
        #                                           overwrite=True)

        # Convolve the VLA data to the common beam:
        com_vla_beam = vla_cube.beams.common_beam()
        print(f"Convolving VLA to common beam: {com_vla_beam}")
        vla_cube = vla_cube.convolve_to(com_vla_beam)
        vla_cube.allow_huge_operations = True

        # Feather with the SD scale factor applied
        feathered_cube = feather_simple_cube(vla_cube.to(u.K),
                                             gbt_cube_specinterp_reproj.to(u.K),
                                             allow_lo_reproj=False,
                                             allow_spectral_resample=False,
                                             allow_huge_operations=True,
                                             lowresscalefactor=this_scfactor,
                                             use_dask=True,
                                             use_save_to_tmp_dir=True)

        # NaN out blank areas post-FFT.
        feathered_cube = feathered_cube.with_mask(vla_cube.mask)

        # this_feathered_filename = vla_data_path / f"{this_vla_filename[:-5]}_feathered.fits"
        this_feathered_filename = vla_data_path / f"{this_vla_filename.name[:-5]}_feathered.K.fits"

        feathered_cube.write(this_feathered_filename, overwrite=True)

        # Divide by the PB.
        pb_cube = SpectralCube.read(vla_data_path / this_vla_filename.name.replace(".fits", "_pb.fits"))

        feathered_cube.allow_huge_operations = True
        pb_cube.allow_huge_operations = True

        feathered_cube_pbcorr = feathered_cube / pb_cube.unitless_filled_data[:]
        this_feathered_filename = vla_data_path / f"{this_vla_filename.name[:-5]}_feathered.K.fits"
        feathered_cube_pbcorr.minimal_subcube().write(this_feathered_filename, overwrite=True)


from glob import glob

for this_linename in hi_cube_name_keys:

    all_filenames = glob(f"*_{this_linename}_*feathered.K.fits")

    for this_name in all_filenames:
        print(this_name)
        # cube = SpectralCube.read(this_name)
        # cube.allow_huge_operations = True
        # pb_cube = SpectralCube.read(this_name.replace("feathered.K.fits", "pb.fits"))
        # pb_cube.allow_huge_operations = True
        # cube_pbcorr = (cube / pb_cube.unitless_filled_data[:]).minimal_subcube()
        this_target, this_config = this_name.split("_")[:2]
        name_stem = this_name.split("_feathered.K.fits")[0]
        output_filename = f"../line_imaging/postprocess/{this_target}/{this_target}_{this_config}+tp_{this_linename}_pbcorr_trimmed_k.fits"
        cube_pbcorr.write(output_filename, overwrite=True)
