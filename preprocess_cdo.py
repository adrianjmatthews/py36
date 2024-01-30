"""Pre preprocess data using cdo.

Before running this script on ada, type 'module add cdo'.
"""

import datetime
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

FILE_MASK=False # Default value

VAR_NAME='swsal'; LEVEL=1062.44; SOURCE='glorys12v1aeq1_zlev_d'

YEAR=2003
#YEAR=range(1998,2022+1)

#MONTH=-999 # if outfile_frequency is 'year'
#MONTH=range(1,12+1) # If outfile_frequency is less than 'year' 
MONTH=1

#------------------------------------------------------------------

iter_year=da.iter_generator(YEAR)
iter_month=da.iter_generator(MONTH)
for year in iter_year:
    for month in iter_month:
        print('### year={0!s} month={1!s}'.format(year,month))
        DIR1=os.path.join(BASEDIR,SOURCE,'raw','tmp')
        DIR2=os.path.join(BASEDIR,SOURCE,'raw')
        xx=VAR_NAME+'_'+str(LEVEL)+'_'+str(year)+str(month).zfill(2)+'.nc'
        file1=os.path.join(DIR1,xx)
        file2=os.path.join(DIR2,xx)
        command='cdo delete,param="irrelevant" '+file1+' '+file2
        print(command)
        os.system(command)
