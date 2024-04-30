"""Average data to lower time resolution using data_analysis.ModifySource."""

import datetime
import os

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt

import data_analysis as da

BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')
#BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

ARCHIVE=True
BASEDIR_ARCHIVE=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

SUBDIR='std'

#VAR_NAME='div'; LEVEL=850; SOURCE1='ncepdoe_plev_6h'; SOURCE2='ncepdoe_plev_d'
#VAR_NAME='uwnd'; LEVEL=850; SOURCE1='erainterim_plev_6h'; SOURCE2='erainterim_plev_d'
VAR_NAME='uwnd'; LEVEL=950; SOURCE1='era5gloerai_plev_h'; SOURCE2='era5gloerai_plev_3h'
#VAR_NAME='ppt'; LEVEL=1; SOURCE1='trmm3b42v7p1_sfc_3h'; SOURCE2='trmm3b42v7p1_sfc_d'
#VAR_NAME='ppt'; LEVEL=1; SOURCE1='imergv07btrm_sfc_30m'; SOURCE2='imergv07btrm_sfc_3h'
#VAR_NAME='ppt'; LEVEL=1; SOURCE1='imergv07btrmrgp1_sfc_3h'; SOURCE2='imergv07btrmrgp1_sfc_d'

#YEAR=2023
#YEAR=range(2001,2022+1)
#MONTH=[-999] # Dummy value if outfile_frequency is 'year'
#MONTH=7
MONTH=range(1,12+1) # If outfile_frequency is less than 'year' 

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

iter_year=da.iter_generator(YEAR)
iter_month=da.iter_generator(MONTH)
for year in iter_year:
    for month in iter_month:
        print('### year={0!s} month={1!s}'.format(year,month))
        aa.year=year
        aa.month=month
        aa.f_time_average(method=3)

if PLOT:
    if True:
        tcoord=aa.cube_out.coord('time')[0]
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
