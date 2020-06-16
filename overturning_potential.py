"""Calculate overturning potential using data_analysis.CubeDiagnostics."""

import datetime
import os

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import pdb

import data_analysis as da
import info

#BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')
BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

ARCHIVE=False
BASEDIR_ARCHIVE=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

LEVEL=500; 
SOURCE='erainterim_plev_d'

FILEPRE='' # e.g., '', '_rac', '_rac_b20_200_n241', '_rac_rm5_n5'
YEAR_BEG=1998; YEAR_END=1998
MONTH1=MONTH2=-999 # Set both MONTH1 and MONTH2 to same (irrelevant) value if outfile_frequency is 'year'
#MONTH1=1; MONTH2=1 # Set month ranges if outfile_frequency is less than 'year'

PLOT=True

VERBOSE=2

#------------------------------------------------------------------

descriptor={}
descriptor['verbose']=VERBOSE
descriptor['basedir']=BASEDIR
descriptor['archive']=ARCHIVE
descriptor['basedir_archive']=BASEDIR_ARCHIVE
descriptor['level']=LEVEL
descriptor['source']=SOURCE
descriptor['filepre']=FILEPRE

# Create instance of CubeDiagnostics object
aa=da.CubeDiagnostics(**descriptor)

# Read psioc data
aa.f_read_data('omega',aa.level)

for year in range(YEAR_BEG,YEAR_END+1):
    for month in range(MONTH1,MONTH2+1):
        print('### year={0!s} month={1!s}'.format(year,month))
        aa.year=year
        aa.month=month
        aa.f_overturning_potential()

if PLOT:
    fig=plt.figure()
    x1=aa.mu
    time_coord=x1.coord('time')
    timec=time_coord.units.num2date(time_coord.points[0])
    x2=x1.extract(iris.Constraint(time=timec))
    #x2=x1
    qplt.contourf(x2)
    plt.gca().coastlines()
    
    plt.show()
    fig.savefig('/gpfs/home/e058/tmp/fig1.png')
