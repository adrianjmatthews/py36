"""Interpolate data to higher time resolution using data_analysis.ModifySource."""

import os

import cftime
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt

import data_analysis as da

#BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')
BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

ARCHIVE=False
BASEDIR_ARCHIVE=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

SUBDIR='std'

VAR_NAME='sst'; LEVEL=1; SOURCE1='sstrey_sfc_7d'; SOURCE2='sstrey_sfc_d'

YEAR_BEG=1981; YEAR_END=2018

PLOT=False

VERBOSE=2

#------------------------------------------------------------------

descriptor={}
descriptor['verbose']=VERBOSE
descriptor['basedir']=BASEDIR
descriptor['archive']=ARCHIVE
descriptor['basedir_archive']=BASEDIR_ARCHIVE
descriptor['subdir']=SUBDIR
descriptor['var_name']=VAR_NAME
descriptor['level']=LEVEL
descriptor['source1']=SOURCE1
descriptor['source2']=SOURCE2

# Create instance of ModifySource object
aa=da.ModifySource(**descriptor)

for year in range(YEAR_BEG,YEAR_END+1):
    print('### year={0!s}'.format(year))
    aa.year=year
    if aa.frequency=='d':
        aa.time1_out=cftime.DatetimeGregorian(aa.year,1,1,0,0)
        aa.time2_out=cftime.DatetimeGregorian(aa.year,12,31,23,59)
    else:
        raise da.ToDoError('Need code for interpolating to other than daily data.')
    aa.f_interpolate_time()

if PLOT:
    #timecon=iris.Constraint(time=cftime.DatetimeGregorian(1985,3,17))
    timecon=iris.Constraint(time=lambda cell: cftime.DatetimeGregorian(1985,1,1)<=cell<=cftime.DatetimeGregorian(1985,12,31))
    latcon=iris.Constraint(latitude=13.5)
    loncon=iris.Constraint(longitude=89.5)
    #x1=aa.cube_out.extract(timecon)
    #
    x1=aa.cube_out.extract(timecon & latcon & loncon)
    x2=aa.cube_in.extract(timecon & latcon & loncon)

    
    #qplt.contourf(x1)
    #plt.gca().coastlines()
    #plt.colorbar()
    #
    qplt.plot(x1)
    qplt.plot(x2)
    
    plt.show()
