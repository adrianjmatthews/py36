"""Regrid data using data_analysis.ModifySource."""

import os

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import pdb

import data_analysis as da

BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')
#BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

ARCHIVE=False
BASEDIR_ARCHIVE=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

# Data to be regridded
#VAR_NAME='ppt'; LEVEL=1; SOURCE1='imergv07btrm_sfc_3h'; SOURCE2='imergv07btrmrg_sfc_3h'
#VAR_NAME='ppt'; LEVEL=1; SOURCE1='trmm3b42v7p1_sfc_d'; SOURCE2='trmm3b42v7p4_sfc_d'
VAR_NAME='uwnd'; LEVEL=950; SOURCE1='era5glo_plev_h'; SOURCE2='era5gloerai_plev_h'
#VAR_NAME='swtheta'; LEVEL=0.494; SOURCE1='glorys12v1eq1_zlev_d'; SOURCE2='glorys12v1eq1erai_zlev_d'

# If SUBDIR is 'std', will regrid data over time
# If SUBDIR is 'processed', will just do a one-off regridding
SUBDIR='std'

# Data set on the target grid. Only the grid is used, not the data
if SOURCE2 in ['notsure1']:
    FILE_GRID=os.path.join(BASEDIR,'ncepdoe_plev_d',SUBDIR,'uwnd_850_2010.nc')
elif SOURCE2 in ['notsure2']:
    FILE_GRID=os.path.join(BASEDIR,SOURCE2,SUBDIR,'dummy_1.nc')
elif SOURCE2 in ['imergv07atrmp1_sfc_3h','imergv07btrmrg_sfc_3h']:
    FILE_GRID=os.path.join(BASEDIR_ARCHIVE,'trmm3b42v7_sfc_3h',SUBDIR,'ppt_1_201912.nc')
elif SOURCE2 in ['era5gloerai_plev_h']:
    FILE_GRID=os.path.join(BASEDIR,'erainterim_plev_6h',SUBDIR,'uwnd_850_1998.nc')
else:
    raise UserWarning('Need to set target grid file for SOURCE2.')

# Option to restrict min/max latitudes on new grid
LATMIN=LATMAX=False # Set both to False to disable this option
#LATMIN=-45; LATMAX=45
#LATMIN=-15; LATMAX=15

#YEAR=2000
#YEAR=range(2001,2003+1)

#MONTH=-999 # if outfile_frequency is 'year'
MONTH=range(1,12+1) # If outfile_frequency is less than 'year' 
#MONTH=7

if SUBDIR=='processed':
    # Set one-off FILEIN1 and FILEOUT1
    dum1='_jjas98-19.nc'
    FILEIN1=os.path.join(BASEDIR,SOURCE1,SUBDIR,VAR_NAME+'_'+str(LEVEL)+dum1)
    FILEOUT1=os.path.join(BASEDIR,SOURCE2,SUBDIR,VAR_NAME+'_'+str(LEVEL)+dum1)

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
descriptor['file_grid']=FILE_GRID
descriptor['latmin']=LATMIN
descriptor['latmax']=LATMAX
if SUBDIR=='processed':
    descriptor['filein1']=FILEIN1
    descriptor['fileout1']=FILEOUT1

# Create instance of ModifySource object
aa=da.ModifySource(**descriptor)

# Get target grid
aa.f_get_target_grid()

iter_year=da.iter_generator(YEAR)
iter_month=da.iter_generator(MONTH)
if SUBDIR=='std':
    for year in iter_year:
        for month in iter_month:
            print('### year={0!s} month={1!s}'.format(year,month))
            aa.year=year
            aa.month=month
            aa.f_regrid()
elif SUBDIR=='processed':
    aa.f_regrid()
else:
    raise UserWarning('Invalid SUBDIR.')

if PLOT:
    fig=plt.figure()

    if True:
        tcoord=aa.cube_in.coord('time')[0]
        time1=tcoord.cell(0)[0]
        timecon=iris.Constraint(time=time1)
        x1=aa.cube_in.extract(timecon)
        qplt.contourf(x1)
        plt.gca().coastlines()
        plt.show()
        fig.savefig('/gpfs/home/e058/tmp/fig1.png')

    if True:
        tcoord=aa.cube_out.coord('time')[0]
        time1=tcoord.cell(0)[0]
        timecon=iris.Constraint(time=time1)
        x2=aa.cube_out.extract(timecon)
        x2=aa.cube_out
        qplt.contourf(x2)
        plt.gca().coastlines()
        plt.show()
        fig.savefig('/gpfs/home/e058/tmp/fig2.png')
