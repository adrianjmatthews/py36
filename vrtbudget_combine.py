"""Combine terms previously calculated from vorticity budget."""

import os

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt

import data_analysis as da

#BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')
BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

ARCHIVE=False
BASEDIR_ARCHIVE=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

LEVEL=850; SOURCE='erainterimEK1_plev_6h'
FILEPRE=''

#YEAR=1998
#YEAR=range(1979,1980+1)
MONTH=[-999] # Dummy value if outfile_frequency is 'year'
#MONTH=range(1,12+1) # If outfile_frequency is less than 'year' 

PLOT=False

VERBOSE=2

#------------------------------------------------------------------

descriptor={}
descriptor['verbose']=VERBOSE
descriptor['basedir']=BASEDIR
descriptor['archive']=ARCHIVE
descriptor['basedir_archive']=BASEDIR_ARCHIVE
descriptor['source']=SOURCE
descriptor['level']=LEVEL
descriptor['filepre']=FILEPRE

# Create instance of CubeDiagnostics object
aa=da.CubeDiagnostics(**descriptor)

# Lazy read data
aa.f_read_data('m_uwnd_dvrtdx',LEVEL)
aa.f_read_data('m_vwnd_dvrtdy',LEVEL)
aa.f_read_data('m_vrt_div',LEVEL)
aa.f_read_data('m_ff_div',LEVEL)
aa.f_read_data('m_domegadx_dvwnddp',LEVEL)
aa.f_read_data('domegady_duwnddp',LEVEL)

iter_year=da.iter_generator(YEAR)
iter_month=da.iter_generator(MONTH)
for year in iter_year:
    for month in iter_month:
        print('### year={0!s} month={1!s}'.format(year,month))
        aa.year=year
        aa.month=month
        aa.f_vrtbudget_combine()

if PLOT:
    tcoord=aa.vrt_stretch.coord('time')
    t1=tcoord.units.num2date(tcoord.points[-1])
    timecon=iris.Constraint(time=t1)
    x1=aa.vrt_stretch.extract(timecon)

    qplt.contourf(x1)
    plt.gca().coastlines()
    
    plt.show()
