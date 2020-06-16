"""Calculate vorticity budget terms."""

import os

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt

import data_analysis as da

BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')

ARCHIVE=False
BASEDIR_ARCHIVE=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

SOURCE='erainterim_plev_6h'
LEVEL_BELOW=550; LEVEL=500; LEVEL_ABOVE=450

#YEAR=1998
MONTH=[-999] # Dummy value if outfile_frequency is 'year'
#MONTH=range(1,12+1) # If outfile_frequency is less than 'year' 

PLOT=False

VERBOSE=2

#------------------------------------------------------------------

descriptor={}
descriptor['verbose']=VERBOSE
descriptor['basedir']=BASEDIR
descriptor['basedir_archive']=BASEDIR_ARCHIVE
descriptor['archive']=ARCHIVE
descriptor['source']=SOURCE
descriptor['level']=LEVEL

# Create instance of CubeDiagnostics object
aa=da.CubeDiagnostics(**descriptor)

# Lazy read data
aa.f_read_data('uwnd',LEVEL_BELOW)
aa.f_read_data('uwnd',LEVEL)
aa.f_read_data('uwnd',LEVEL_ABOVE)
aa.f_read_data('vwnd',LEVEL_BELOW)
aa.f_read_data('vwnd',LEVEL)
aa.f_read_data('vwnd',LEVEL_ABOVE)
aa.f_read_data('vrt',LEVEL_BELOW)
aa.f_read_data('vrt',LEVEL)
aa.f_read_data('vrt',LEVEL_ABOVE)
aa.f_read_data('omega',LEVEL)
aa.f_read_data('div',LEVEL,verbose=True)

iter_year=da.iter_generator(YEAR)
iter_month=da.iter_generator(MONTH)
for year in iter_year:
    for month in iter_month:
        print('### year={0!s} month={1!s}'.format(year,month))
        aa.year=year
        aa.month=month
        aa.f_vrtbudget(LEVEL_BELOW,LEVEL,LEVEL_ABOVE)

if PLOT:
    tcoord=aa.dvrtdt.coord('time')
    time1=tcoord.units.num2date(tcoord.points[-1])
    time_constraint=iris.Constraint(time=time1)
    x1=aa.dvrtdt.extract(time_constraint)

    qplt.contourf(x1)
    plt.gca().coastlines()
    
    plt.show()
