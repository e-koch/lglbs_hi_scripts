
'''
Runs a combination of weighting and array configurations to
return the PSF shape and achieved noise.

Current example uses ic10ctr in HI wtih 1.2 km/s channels.
'''

from pathlib import Path


from casatools import ms
from casatools import imager
from casatools import synthesisutils
from casatools import msmetadata

synthutil = synthesisutils()
myms = ms()

do_config_imaing = False
do_combo_imaging = True


this_field = "ic10ctr"
this_line = "hi21cm_1p2kms"

data_path = Path("/reduction10/erickoch/LGLBS/line_imaging/imaging/ic10ctr/")

output_path = Path("/reduction10/erickoch/LGLBS/line_imaging/psf_and_noise_tests/")

configs = ['A', 'B', 'C', 'D']

config_combos = ['A+B+C+D', 'A+B+C', 'B+C+D', 'C+D']
# config_combos = ['A+B+C']

briggs_weights = [f"briggs_{this_r}" for this_r in
                  [2, 1.5, 1.0, 0.5, 0.0, -0.5, -1.0, -1.5, -2.0]]
weightings = ['natural', 'uniform'] + briggs_weights

all_vis_dict = {}

for this_config in configs:

    this_vis = data_path / f"{this_field}_{this_config}_{this_line}.ms"

    if not this_vis.exists():
        raise FileExistsError(f"Can't find {this_vis.name}")

    all_vis_dict[this_config] = this_vis

# Loop through individual configs:
if do_config_imaing:
    for this_config in configs:

        myvis = all_vis_dict[this_config]

        # Set up imaging:

        this_im = imager()
        this_im.selectvis(vis=str(myvis),
                        field="",
                        spw=str(0))
        image_settings = this_im.advise()
        this_im.close()

        cell_size = [image_settings[2]['value'],
                    image_settings[2]['unit']]

        this_msmd = msmetadata()
        this_msmd.open(str(myvis))
        nchan = this_msmd.nchan(0)
        mean_freq = this_msmd.chanfreqs(int(0)).mean() / 1.e9 # Hz to GHz
        this_msmd.close()

        approx_pbsize = 1.5 * (45. / mean_freq) * 60 # arcsec
        this_imsize = synthutil.getOptimumSize(int(approx_pbsize / image_settings[2]['value']))


        this_cellsize = f"{round(cell_size[0] * 0.8, 1)}{cell_size[1]}"

        this_pblim = 0.25

        this_nsigma = 5.
        this_niter = 0

        # 0.05 of the total range. The idea is to find a signal free channel.
        this_start = int(0.05 * nchan)

        for this_weight in weightings:

            this_imagename = f"{this_field}_{this_config}_{this_line}_{this_weight}"

            # Clean up any possible imaging remnants first
            rmtables(f"{this_imagename}*")

            tclean(vis=str(myvis),
                    field="",
                    spw=str(0),
                    cell=this_cellsize,
                    imsize=this_imsize,
                    specmode='cube',
                    weighting=this_weight.split("_")[0],
                    robust=float(this_weight.split("_")[1]) if "briggs" in this_weight else 0.0,
                    start=this_start,
                    width=1,
                    nchan=1,
                    niter=this_niter,
                    nsigma=this_nsigma,
                    imagename=str(this_imagename),
                    restfreq="1.42040575177GHz",
                    pblimit=this_pblim)

            # Estimate the expected sensitivity
            # only mfs is supported right now, so skip
            # out = apparentsens(myvis,
            #                     field=this_field,
            #                     spw=str(0),
            #                     cell=this_cellsize,
            #                     imsize=this_imsize,
            #                     specmode='mfs',
            #                     weighting='briggs',
            #                     robust=0.0)


# Now loop through the array combos
if do_combo_imaging:

    for this_combo in config_combos:

        these_vis = [all_vis_dict[vis] for vis in this_combo.split('+')]


        myvis = all_vis_dict[this_combo.split("+")[0]]

        # Set up imaging:

        this_im = imager()
        this_im.selectvis(vis=str(myvis),
                        field="",
                        spw=str(0))
        image_settings = this_im.advise()
        this_im.close()

        cell_size = [image_settings[2]['value'],
                    image_settings[2]['unit']]

        this_msmd = msmetadata()
        this_msmd.open(str(myvis))
        nchan = this_msmd.nchan(0)
        mean_freq = this_msmd.chanfreqs(int(0)).mean() / 1.e9 # Hz to GHz
        this_msmd.close()

        approx_pbsize = 1.5 * (45. / mean_freq) * 60 # arcsec
        this_imsize = synthutil.getOptimumSize(int(approx_pbsize / image_settings[2]['value']))


        this_cellsize = f"{round(cell_size[0] * 0.8, 1)}{cell_size[1]}"

        this_pblim = 0.25

        this_nsigma = 5.
        this_niter = 0

        # 0.05 of the total range. The idea is to find a signal free channel.
        this_start = int(0.05 * nchan)

        for this_weight in weightings:

            this_imagename = f"{this_field}_{this_combo}_{this_line}_{this_weight}"

            # Clean up any possible imaging remnants first
            rmtables(f"{this_imagename}*")

            tclean(vis=[str(myvis) for myvis in these_vis],
                    field="",
                    spw=str(0),
                    cell=this_cellsize,
                    imsize=this_imsize,
                    specmode='cube',
                    weighting=this_weight.split("_")[0],
                    robust=float(this_weight.split("_")[1]) if "briggs" in this_weight else 0.0,
                    start=this_start,
                    width=1,
                    nchan=1,
                    niter=this_niter,
                    nsigma=this_nsigma,
                    imagename=str(this_imagename),
                    restfreq="1.42040575177GHz",
                    pblimit=this_pblim)
