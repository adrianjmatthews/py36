"""Calculate psi, chi, vrt, div, wndspeed from uwnd, vwnd using data_analysis.Wind."""

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

SUBDIR='std' # 'std': input data with time axis
                   # 'processed': input data with no time axis

# Optional
#TDOMAINID='djf7980-1718'

#LEVEL=1000; 
SOURCE='erainterim_plev_d'
#LEVEL=200; SOURCE='ncepdoe_plev_d'
#LEVEL=50; SOURCE='ncepncar_plev_d'
#LEVEL=10; SOURCE='ncepdoe_zlev_d'

FLAG_PSI=False
FLAG_CHI=False
FLAG_VRT=False
FLAG_DIV=False
FLAG_WNDSPD=False
FLAG_UWNDCHI=True
FLAG_VWNDCHI=True

if SUBDIR=='std':
    FILEPRE='' # e.g., '', '_rac', '_rac_b20_200_n241', '_rac_rm5_n5'
    YEAR_BEG=1999; YEAR_END=1999
    MONTH1=MONTH2=-999 # Set both MONTH1 and MONTH2 to same (irrelevant) value if outfile_frequency is 'year'
    #MONTH1=1; MONTH2=1 # Set month ranges if outfile_frequency is less than 'year'
elif SUBDIR=='processed':
    FILEPRE='_'+TDOMAINID # e.g., TDOMAINID of time mean data.
else:
    raise ValueError('SUBDIR is invalid.')

PLOT=False

VERBOSE=2

#------------------------------------------------------------------

descriptor={}
descriptor['verbose']=VERBOSE
descriptor['basedir']=BASEDIR
descriptor['archive']=ARCHIVE
descriptor['basedir_archive']=BASEDIR_ARCHIVE
descriptor['subdir']=SUBDIR
descriptor['source']=SOURCE
descriptor['level']=LEVEL
descriptor['filepre']=FILEPRE
descriptor['flag_psi']=FLAG_PSI
descriptor['flag_chi']=FLAG_CHI
descriptor['flag_vrt']=FLAG_VRT
descriptor['flag_div']=FLAG_DIV
descriptor['flag_wndspd']=FLAG_WNDSPD
descriptor['flag_uwndchi']=FLAG_UWNDCHI
descriptor['flag_vwndchi']=FLAG_VWNDCHI

# Create instance of Wind object
aa=da.Wind(**descriptor)

if SUBDIR=='std':
    for year in range(YEAR_BEG,YEAR_END+1):
        for month in range(MONTH1,MONTH2+1):
            print('### year={0!s} month={1!s}'.format(year,month))
            aa.year=year
            aa.month=month
            aa.f_wind()
elif SUBDIR=='processed':
    aa.f_wind()

if PLOT:
    x1=aa.vwndchi
    time_coord=x1.coord('time')
    timec=time_coord.units.num2date(time_coord.points[0])
    x2=x1.extract(iris.Constraint(time=timec))
    #x2=x1
    qplt.contourf(x2)
    plt.gca().coastlines()
    
    plt.show()
