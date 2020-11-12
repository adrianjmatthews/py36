"""Preprocess data using data_analysis.DataConverter."""

import datetime
import os

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt

import data_analysis as da

BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')
#BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')

ARCHIVE=False
BASEDIR_ARCHIVE=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

FILE_MASK=False # Default value

VAR_NAME='uwnd'; LEVEL=925; SOURCE='era5trp_plev_h'
#VAR_NAME='psfc'; LEVEL=1; SOURCE='erainterim_sfc_d'
#VAR_NAME='vrt'; LEVEL=500; SOURCE='erainterim_plev_6h'
#VAR_NAME='zg'; LEVEL=250; SOURCE='hadgem2esajhog_plev_d'
#VAR_NAME='uwnd'; LEVEL=1000; SOURCE='ncepdoe_plev_d'
#VAR_NAME='psfc'; LEVEL=1; SOURCE='ncepdoe_sfc_d'
#VAR_NAME='vwnd'; LEVEL=10; SOURCE='ncepdoegg_zlev_d'
#VAR_NAME='vwnd'; LEVEL=1; SOURCE='ncepncar_sfc_d'
#VAR_NAME='ta'; LEVEL=200; SOURCE='ncepncar_plev_d'
#VAR_NAME='olr'; LEVEL=0; SOURCE='olrcdr_toa_d'
#VAR_NAME='olr'; LEVEL=0; SOURCE='olrinterp_toa_d'
#VAR_NAME='tsc'; LEVEL='all'; SOURCE='sg579m031oi01_zlev_h'
#VAR_NAME='sst'; LEVEL=1; SOURCE='sstrey_sfc_7d'; FILE_MASK='lsmask.nc'
#VAR_NAME='tauy'; LEVEL=1; SOURCE='tropflux_sfc_d'
#VAR_NAME='ppt'; LEVEL=1; SOURCE='trmm3b42v7_sfc_3h'
#VAR_NAME='ppt'; LEVEL=1; SOURCE='imergmt2_sfc_30m'

YEAR_BEG=2018; YEAR_END=2018 # if outfile_frequency is 'year' or less

#MONTH1=MONTH2=-999 # if outfile_frequency is 'year'
MONTH1=11; MONTH2=11 # if outfile_frequency is less than 'year'

PLOT=True

VERBOSE=2

#------------------------------------------------------------------

descriptor={}
descriptor['verbose']=VERBOSE
descriptor['basedir']=BASEDIR
descriptor['archive']=ARCHIVE
descriptor['basedir_archive']=BASEDIR_ARCHIVE
descriptor['var_name']=VAR_NAME
descriptor['level']=LEVEL
descriptor['source']=SOURCE
descriptor['file_mask']=FILE_MASK

aa=da.DataConverter(**descriptor)

for year in range(YEAR_BEG,YEAR_END+1):
    for month in range(MONTH1,MONTH2+1):
        print('### year={0!s} month={1!s}'.format(year,month))
        aa.year=year
        aa.month=month
        aa.read_cube()
        aa.format_cube()
        aa.write_cube()

if PLOT:
    tcoord=aa.cube.coord('time')
    time1=tcoord.units.num2date(tcoord.points[0])
    time_constraint=iris.Constraint(time=time1)
    x1=aa.cube.extract(time_constraint)

    #x1=aa.cube
    qplt.contourf(x1)
    plt.gca().coastlines()
    
    plt.show()
