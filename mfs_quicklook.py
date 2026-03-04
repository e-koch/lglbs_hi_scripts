os.system('rm -rf '+outfile_root+'.*')

tclean(vis=infile,
       imagename=outfile_root,
       phasecenter=phase_center,
       imsize=imsize,
       cell=cell,
       datacolumn='data',
       specmode='mfs',
       nterms=1,
       gridder='mosaic',
       niter=0,
       weighting='briggs',
       robust=0.5,
       )
       
