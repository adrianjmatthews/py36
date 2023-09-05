"""Preprocess data using data_analysis.DataConverter."""

import datetime
import os

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt

import data_analysis as da

#BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data') # UEA
BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data') # UEA
#BASEDIR=os.path.join(os.path.sep,'gws','nopw','j04','bobble','matthews','data') # JASMIN

ARCHIVE=False
BASEDIR_ARCHIVE=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

FILE_MASK=False # Default value

#VAR_NAME='ppt'; LEVEL=1; SOURCE='cmap_sfc_5d'
#VAR_NAME='uwnd'; LEVEL=975; SOURCE='era5mcw_plev_h'
#VAR_NAME='vrt'; LEVEL=825; SOURCE='era5glo_plev_h'
#VAR_NAME='psfc'; LEVEL=1; SOURCE='erainterim_sfc_d'
#VAR_NAME='zg'; LEVEL=850; SOURCE='erainterim_plev_6h'
#VAR_NAME='zg'; LEVEL=250; SOURCE='hadgem2esajhog_plev_d'
#VAR_NAME='ppt'; LEVEL=1; SOURCE='imergtrm_sfc_30m'
#VAR_NAME='zg'; LEVEL=200; SOURCE='ncepdoe_plev_d'
#VAR_NAME='psfc'; LEVEL=1; SOURCE='ncepdoe_sfc_d'
#VAR_NAME='vwnd'; LEVEL=10; SOURCE='ncepdoegg_zlev_d'
#VAR_NAME='vwnd'; LEVEL=1; SOURCE='ncepncar_sfc_d'
#VAR_NAME='ta'; LEVEL=200; SOURCE='ncepncar_plev_d'
#VAR_NAME='olr'; LEVEL=0; SOURCE='olrcdr_toa_d'
#VAR_NAME='olr'; LEVEL=0; SOURCE='olrinterp_toa_d'
#VAR_NAME='ssft'; LEVEL=1; SOURCE='ostial4reptrp_sfc_d'
#VAR_NAME='tsc'; LEVEL='all'; SOURCE='sg579m031oi01_zlev_h'
#VAR_NAME='sst'; LEVEL=1; SOURCE='sstrey_sfc_7d'; FILE_MASK='lsmask.nc'
#VAR_NAME='tauy'; LEVEL=1; SOURCE='tropflux_sfc_d'
#VAR_NAME='ppt'; LEVEL=1; SOURCE='trmm3b42v7_sfc_3h'
VAR_NAME='swtheta'; LEVEL=0.494; SOURCE='glorys12v1eq1_zlev_d'

YEAR=2017
#YEAR=range(1998,2022+1)

#MONTH=-999 # if outfile_frequency is 'year'
#MONTH=range(1,12+1) # If outfile_frequency is less than 'year' 
MONTH=1

PLOT=False

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

iter_year=da.iter_generator(YEAR)
iter_month=da.iter_generator(MONTH)
for year in iter_year:
    for month in iter_month:
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
