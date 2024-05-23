"""Calculate surface wind stress component from wind component."""

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

LEVEL=1; SOURCE='era5gloerai_sfc_d'

METHOD='uwnd' # 'uwnd', 'vwnd', or 'both'

YEAR=1998
#YEAR=range(1998,2023+1)

#MONTH=-999 # if outfile_frequency is 'year'
MONTH=1
#MONTH=range(1,12+1) # If outfile_frequency is less than 'year' 

PLOT=False

VERBOSE=2
#------------------------------------------------------------------

descriptor={}
descriptor['verbose']=VERBOSE
descriptor['basedir']=BASEDIR
descriptor['source']=SOURCE
descriptor['level']=LEVEL
descriptor['archive']=ARCHIVE

# Create instance of CubeDiagnostics object
aa=da.CubeDiagnostics(**descriptor)

# Lazy read data
if METHOD not in ['uwnd','vwnd','both']:
    raise UserWarning('Invalid METHOD.')
if METHOD in['uwnd','both']:
    aa.f_read_data('uwnd',LEVEL)
if METHOD in['vwnd','both']:
    aa.f_read_data('vwnd',LEVEL)

iter_year=da.iter_generator(YEAR)
iter_month=da.iter_generator(MONTH)
for year in iter_year:
    for month in iter_month:
        print('### year={0!s} month={1!s}'.format(year,month))
        aa.year=year
        aa.month=month
        aa.f_tau_from_wind(METHOD)

if PLOT:
    fig=plt.figure()

    tcoord=aa.sa.coord('time')
    t1=tcoord.units.num2date(tcoord.points[-1])
    timecon=iris.Constraint(time=t1)
    x1=aa.sa.extract(timecon)

    qplt.contourf(x1)
    plt.gca().coastlines()
    
    plt.show()

    fig.savefig('/gpfs/home/e058/tmp/fig1.png')
