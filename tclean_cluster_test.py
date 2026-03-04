
'''
Writing a cluster-friendly version of tclean that allows
for estimating when to exit for simple restarting of jobs.


'''

# import os
# import shutil
# import copy
import time

# import numpy

from casatasks.private.imagerhelpers.imager_base import PySynthesisImager
from casatasks.private.imagerhelpers.imager_parallel_continuum import PyParallelContSynthesisImager
from casatasks.private.imagerhelpers.imager_parallel_cube import PyParallelCubeSynthesisImager

from casatasks.private.imagerhelpers.input_parameters import ImagerParameters



def tclean_restartable(vis, imagename, *args, **kwargs):

    # (2) Set up Input Parameters
    # - List all parameters that you need here
    # - Defaults will be assumed for unspecified parameters
    # - Nearly all parameters are identical to that in the task. Please look at the
    # list of parameters under __init__ using " help ImagerParameters " )

    msname = vis

    kwargs['loopgain'] = kwargs.pop('gain', 0.1)
    kwargs['timestr'] = kwargs.pop('timerange', '')
    kwargs['dopbcorr'] = kwargs.pop('pbcor', False)
    kwargs['uvdist'] = kwargs.pop('uvrange', '')
    kwargs['scalebias'] = kwargs.pop('smallscalebias', 0.0)
    kwargs['cfcache'] = kwargs.pop('cfcache', '')
    kwargs['nchan'] = kwargs.pop('nchan', -1) or -1
    kwargs['nterms'] = kwargs.pop('nterms', 2) or 2
    kwargs['facets'] = kwargs.pop('facets', 1) or 1
    kwargs['chanchunks'] = kwargs.pop('chanchunks', -1) or -1

    kwargs.pop('calcpsf', True)
    kwargs.pop('psfphasecenter', None)
    kwargs.pop('calcres', True)
    kwargs.pop('restoration', True)
    kwargs.pop('intent', "")
    kwargs.pop('observation', "")
    kwargs.pop('selectdata', True)

    for key in kwargs.keys():
        if kwargs[key] is None:
            kwargs.pop(key)

    paramList = ImagerParameters(msname=msname, imagename=imagename, **kwargs)

    # (3) Construct the PySynthesisImager object, with all input parameters

    imager = PySynthesisImager(params=paramList)

    # (4) Initialize various modules.
    # - Pick only the modules you will need later on. For example, to only make
    # the PSF, there is no need for the deconvolver or iteration control modules.

    # Initialize modules major cycle modules

    imager.initializeImagers()
    imager.initializeNormalizers()
    imager.setWeighting()

    # Init minor cycle modules

    imager.initializeDeconvolvers()
    imager.initializeIterationControl()

    # (5) Make the initial images

    imager.makePSF()
    imager.makePB()
    imager.runMajorCycle()  # Make initial dirty / residual image

    # (6) Make the initial clean mask
    imager.hasConverged()
    imager.updateMask()

    majornumber = 0
    minornumber = 0

    majornumber_max = 3

    # (7) Run the iteration loops
    # while (not imager.hasConverged()) or (majornumber > majornumber_max):
    while (majornumber <= majornumber_max):

        casalog.post(f"***Niter for this cycle is: {imager.iterpars['niter']}" , "INFO", "task_tclean");

        print(imager.iterpars['niter'])
        print(imager.IBtool.getiterationdetails())

        imager.runMinorCycle()
        # exportfits(imagename+".image.tt0", imagename+".image.tt0.minor{0}.fits".format(minornumber), overwrite=True)
        # exportfits(imagename+".residual.tt0", imagename+".residual.tt0.minor{0}.fits".format(minornumber), overwrite=True)
        minornumber += 1
        imager.runMajorCycle()
        # exportfits(imagename+".image.tt0", imagename+".image.tt0.major{0}.fits".format(majornumber), overwrite=True)
        # exportfits(imagename+".residual.tt0", imagename+".residual.tt0.minor{0}.fits".format(minornumber), overwrite=True)
        # exportfits(imagename+".model.tt0", imagename+".model.tt0.minor{0}.fits".format(minornumber), overwrite=True)
        majornumber += 1
        imager.updateMask()

        # Update the niter per cycle
        imager.iterpars['niter'] *= 5

        print(imager.IBtool.cleanComplete())

    # (8) Finish up

    retrec = imager.getSummary()
    imager.restoreImages()

    # imager.pbcorImages()

    # (9) Close tools.

    imager.deleteTools()

    return retrec



####
# Testing example
####

myvis = "orig_ms_files/Brick1Tile1_2019.1.01182.S_uid___A002_Xe29133_X1e6a_12CO21.ms.contsub"

rmtables(f'{myvis}.tclean_nitervariation_test*')

retrec = tclean_restartable(vis=myvis,
       imagename=f'{myvis}.tclean_nitervariation_test',
       specmode='cube', cell='1arcsec', imsize=1024,
       niter=1000, nsigma=2.,
       start='-300km/s', width='5km/s',nchan=40,
       weighting='briggs', robust=0.5,
       usemask='pb', pbmask=0.3,
       pblimit=0.3)
