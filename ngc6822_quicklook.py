inroot = '../lglbs/imaging/ngc6822/'

phase_center = 'J2000 19h44m57.74s -14d48m12.4s'
cell = '5arcsec'
imsize = [1000,1000]

infile = inroot+'ngc6822_20A-346_C_1_hi21cm.ms.contsub'
outfile_root = inroot+'ngc6822_C_1_quicklook'
execfile('mfs_quicklook.py')

infile = inroot+'ngc6822_20A-346_C_2_hi21cm.ms.contsub'
outfile_root = inroot+'ngc6822_C_2_quicklook'
execfile('mfs_quicklook.py')

infile = inroot+'ngc6822_20A-346_C_3_hi21cm.ms.contsub'
outfile_root = inroot+'ngc6822_C_3_quicklook'
execfile('mfs_quicklook.py')

infile = inroot+'ngc6822_20A-346_C_4_hi21cm.ms.contsub'
outfile_root = inroot+'ngc6822_C_4_quicklook'
execfile('mfs_quicklook.py')
