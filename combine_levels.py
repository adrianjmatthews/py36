"""Combine levels into single cube using data_analysis.CombineLevels"""

import datetime
import os

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt

import data_analysis as da
import info

#BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')
BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

ARCHIVE=False
BASEDIR_ARCHIVE=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

SUBDIR='processed' # 'std': input data with time axis
                   # 'processed': input data with no time axis

# Optional
TDOMAINID='ann1999'

VAR_NAME='omega_lambda';SOURCE='ncepdoe_plev_d'
LEVELS=info.ncepdoe_levels

if SUBDIR=='std':
    FILEPRE='' # e.g., '', '_rac', '_rac_b20_200_n241', '_rac_rm5_n5'
    YEAR_BEG=1999; YEAR_END=1999
    MONTH1=MONTH2=-999 # Set both MONTH1 and MONTH2 to same (irrelevant) value if outfile_frequency is 'year'
    #MONTH1=1; MONTH2=1 # Set month ranges if outfile_frequency is less than 'year'
elif SUBDIR=='processed':
    FILEPRE='_'+TDOMAINID # e.g., TDOMAINID of time mean data.
else:
    raise ValueError('SUBDIR is invalid.')

PLOT=True

VERBOSE=2

#------------------------------------------------------------------

descriptor={}
descriptor['verbose']=VERBOSE
descriptor['basedir']=BASEDIR
descriptor['archive']=ARCHIVE
descriptor['basedir_archive']=BASEDIR_ARCHIVE
descriptor['subdir']=SUBDIR
descriptor['source']=SOURCE
descriptor['var_name']=VAR_NAME
descriptor['levels']=LEVELS
descriptor['filepre']=FILEPRE


# Create instance of CombineLevels object
aa=da.CombineLevels(**descriptor)

if SUBDIR=='std':
    for year in range(YEAR_BEG,YEAR_END+1):
        for month in range(MONTH1,MONTH2+1):
            print('### year={0!s} month={1!s}'.format(year,month))
            aa.year=year
            aa.month=month
            aa.f_combine_levels()
elif SUBDIR=='processed':
    aa.f_combine_levels()

if PLOT:
    fig=plt.figure()
    x1=aa.data_all
    x2=x1.extract(iris.Constraint(latitude=0.0))
    qplt.contourf(x2,yrev=1)
    
    plt.show()
    fig.savefig('/gpfs/home/e058/tmp/fig1.png')
