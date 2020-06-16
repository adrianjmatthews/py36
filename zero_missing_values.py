"""Set missing values to zero using data_analysis.ModifySource."""

import datetime
import os

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt

import data_analysis as da

#BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')
BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

ARCHIVE=False
BASEDIR_ARCHIVE=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

SUBDIR='std'

VAR_NAME='ppt'; LEVEL=1; SOURCE1='trmm3b42v7_sfc_3h'; SOURCE2='trmm3b42v7p1_sfc_3h'

YEAR_BEG=2019; YEAR_END=2019 # if outfile_frequency is 'year' or less

#MONTH1=MONTH2=-999 # if outfile_frequency is 'year'
MONTH1=5; MONTH2=9 # if outfile_frequency is less than 'year'

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

# Create instance of ModifySource object
aa=da.ModifySource(**descriptor)

for year in range(YEAR_BEG,YEAR_END+1):
    for month in range(MONTH1,MONTH2+1):
        print('### year={0!s} month={1!s}'.format(year,month))
        aa.year=year
        aa.month=month
        aa.f_zero_missing_values()

if PLOT:
    if True:
        tcoord=aa.cube_out.coord('time')[0]
        #time1=tcoord.units.num2date(tcoord.cell(0)[0])
        time1=tcoord.cell(0)[0]
        timecon=iris.Constraint(time=time1)
        x1=aa.cube_out.extract(timecon)
        qplt.contourf(x1)
        plt.gca().coastlines()
    elif False:
        #latcon=iris.Constraint(latitude=0.125)
        #loncon=iris.Constraint(longitude=80.125)
        latcon=iris.Constraint(latitude=50.0)
        loncon=iris.Constraint(longitude=140.0)
        x1=aa.cube_in.extract(latcon & loncon)
        x2=aa.cube_out.extract(latcon & loncon)
        qplt.plot(x1,label='in')
        qplt.plot(x2,label='out')
        plt.legend()
    
    plt.show()
