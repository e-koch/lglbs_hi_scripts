
import Base.Threads.@spawn

using CloudClean, FITSIO
using ProgressBars

<<<<<<< HEAD
f = FITS("NGC6822-center_cube_ircs.fits")
=======
cube_name, mask_name, finite_mask_name = parse.(String, ARGS);

# cube_name = "NGC6822-center_cube_ircs.fits"
# mask_name = "NGC6822-center_cube_ircs_nosourcemask.fits"
# finite_mask_name = "NGC6822-center_cube_ircs_finitemask.fits"

f = FITS(cube_name)
>>>>>>> 24ab4b9 (Long overdue update, including foreground separation, feathering scripts, rms estimation, and cluster imaging)
out_image = read(f[1]);
header = read_header(f[1]);
close(f)


<<<<<<< HEAD
f = FITS("NGC6822-center_cube_ircs_nosourcemask.fits")
=======
f = FITS(mask_name)
>>>>>>> 24ab4b9 (Long overdue update, including foreground separation, feathering scripts, rms estimation, and cluster imaging)
out_mask = read(f[1]);
close(f)

# This mask only includes where the GBT map covers
<<<<<<< HEAD
f = FITS("NGC6822-center_cube_ircs_finitemask.fits")
=======
f = FITS(finite_mask_name)
>>>>>>> 24ab4b9 (Long overdue update, including foreground separation, feathering scripts, rms estimation, and cluster imaging)
out_mask_finite = read(f[1]);
close(f)

galactic_model_cube = similar(out_image, Float32);

# Define parameters for cloudclean
# NOTE: these may not be optimal for your data!
Np = 51;
widx = 129;
seed = 2021;
# ndraw = 10

Threads.@threads for k in ProgressBar(1:size(out_image, 3))

    this_channel = out_image[:, :, k];
    this_channel_mask = Bool.(out_mask[:, :, k]);
    this_channel_mask_finite = Bool.(out_mask_finite[:, :, k]);

    # If there's no masking within the GBT footprint, skip
    if this_channel_mask == this_channel_mask_finite
        println("Skipping channel $k")
        galactic_model_cube[:, :, k] = this_channel;
        continue
    end

    result = proc_continuous(this_channel,
                            .!this_channel_mask,
                            Np=Np,
                            widx=widx,
                            seed=seed,);
                            # ndraw=10)

    # Record the meaned model
    galactic_model_cube[:, :, k] = result;

end

f = FITS("NGC6822-center_cube_ircs_galactic_model_cloudclean_noheader.fits", "w")
write(f, galactic_model_cube);
close(f)



