
'''
This script is for identifying problematic scans that persist across

'''
fields = ['IC1613_1_CTR',]
	 # 'IC1613_2', 'IC1613_3',
         # 'IC1613_4', 'IC1613_5', 'IC1613_6',
         # 'IC1613_7']

myvis = '20A-346.sb39983867.eb40024999.59433.251063668984.speclines.ms.split'

#myvis = '20A-346.sb40427463.eb40926917.59518.03342894676.speclines.ms.split'
#myvis = '20A-346.sb40427463.eb40928530.59519.01410726852.speclines.ms.split'

myspw = '2'
# myspw = '5'

# Velocity
#mystart = '-280km/s'
#mynchan = 20
#mywidth ='5km/s'

# Channel
mystart = 2000
mynchan = 100
mywidth = 1

#mycell = '1arcsec'
mycell = '3arcsec'

# Small cube with all scans to check which fields have problematic
# data

for this_field in fields:

    # for this_robust in [-2., -1., 0., 0.5, 1., 2.]:
    for this_robust in [0.5, 1., 2.]:

	    outname = f"{myvis}.{this_field}.spw{myspw}.perchan.robust{this_robust:.1f}.ea11ea20_flag"
	    #outname = f"{myvis}.{this_field}.spw{myspw}.addedflagging"

	    rmtables(f"{outname}.*")

	    tclean(vis=myvis,
        	   field=this_field,
	           imagename=outname,
	           imsize=1024,
	           cell=mycell,
	           datacolumn='data',
	           specmode='cube',
	           gridder='mosaic',
	           spw=myspw,
	           start=mystart,
	           nchan=mynchan,
	           width=mywidth,
        	   niter=0,
	           weighting='briggs',
	           robust=this_robust,
               restfreq='1.420405GHz',
        	   )
