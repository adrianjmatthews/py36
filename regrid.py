"""Regrid data using data_analysis.ModifySource."""

import datetime
import os

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import pdb

import data_analysis as da

#BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')
BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

ARCHIVE=False
BASEDIR_ARCHIVE=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

# Data to be regridded
VAR_NAME='ppt'; LEVEL=1; SOURCE1='trmm3b42v7p1_sfc_d'; SOURCE2='trmm3b42v7p4_sfc_d'

# If SUBDIR is 'std', will regrid data over time
# If SUBDIR is 'processed', will just do a one-off regridding
SUBDIR='processed'

# Data set on the target grid. Only the grid is used, not the data
#FILE_GRID=os.path.join(BASEDIR,'ncepdoe_plev_d',SUBDIR,'uwnd_850_2010.nc')
FILE_GRID=os.path.join(BASEDIR,SOURCE2,SUBDIR,'dummy_1.nc')

# Option to restrict min/max latitudes on new grid
#LATMIN=LATMAX=False # Set both to False to disable this option
LATMIN=-45; LATMAX=45


YEAR_BEG=1998; YEAR_END=2018 # if outfile_frequency is 'year' or less

#MONTH1=MONTH2=-999 # if outfile_frequency is 'year'
MONTH1=1; MONTH2=12 # if outfile_frequency is less than 'year'

if SUBDIR=='processed':
    # Set one-off FILEIN1 and FILEOUT1
    dum1='_jjas98-12.nc'
    FILEIN1=os.path.join(BASEDIR,SOURCE1,SUBDIR,VAR_NAME+'_'+str(LEVEL)+dum1)
    FILEOUT1=os.path.join(BASEDIR,SOURCE2,SUBDIR,VAR_NAME+'_'+str(LEVEL)+dum1)

PLOT=True

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

if SUBDIR=='std':
    for year in range(YEAR_BEG,YEAR_END+1):
        for month in range(MONTH1,MONTH2+1):
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

    if False:
        tcoord=aa.cube_in.coord('time')[0]
        time1=tcoord.cell(0)[0]
        timecon=iris.Constraint(time=time1)
        x1=aa.cube_in.extract(timecon)
        qplt.contourf(x1)
        plt.gca().coastlines()
        plt.show()
        fig.savefig('/gpfs/home/e058/tmp/fig1.png')

    if True:
        #tcoord=aa.cube_out.coord('time')[0]
        #time1=tcoord.cell(0)[0]
        #timecon=iris.Constraint(time=time1)
        #x2=aa.cube_out.extract(timecon)
        x2=aa.cube_out
        qplt.contourf(x2)
        plt.gca().coastlines()
        plt.show()
        fig.savefig('/gpfs/home/e058/tmp/fig2.png')
