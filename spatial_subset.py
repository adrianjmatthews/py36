"""Create a spatially subsetted field using data_analysis.SpatialSubset."""

import datetime
import os

import iris.quickplot as qplt
import matplotlib.pyplot as plt
import pdb

import data_analysis as da

#BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')
BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

ARCHIVE=False
BASEDIR_ARCHIVE=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

SUBDIR='std' # 'std': input data with time axis
                   # 'processed': input data with no time axis

# Optional
#TDOMAINID='CCEK102E98-18-00UTC-and-M0002b'
#TDOMAINID='M0002b-and-not-CCEK102E98-18-00UTC'
#TDOMAINID='CCEK102E01-20-00UTC-and-M0002a'
#TDOMAINID='M0002a-and-not-CCEK102E01-20-00UTC'

#VAR_NAME='div'; LEVEL=975; SOURCE='erainterim_plev_6h'
#VAR_NAME='ta'; LEVEL=1; SOURCE='era5plp_sfc_h'
#VAR_NAME='omega'; LEVEL=400; SOURCE='ncepdoe_plev_d'
#VAR_NAME='vwndptap'; LEVEL=50; SOURCE='ncepncar_plev_d'
#VAR_NAME='olr'; LEVEL=0; SOURCE='olrinterp_toa_d'
#VAR_NAME='sst'; LEVEL=1; SOURCE='sstrey_sfc_d'
#VAR_NAME='ppt'; LEVEL=1; SOURCE='trmm3b42v7p1_sfc_3h'
VAR_NAME='ppt'; LEVEL=1; SOURCE='imergmcw_sfc_30m'

# Compulsory band1
BAND1_NAME='latitude' # Dimension to average over: 'latitude' or 'longitude'
#BAND1_VAL1=-2.625 # TRMM precip for CCKW trajectories
#BAND1_VAL1=-15 # TRMM precip for CCER trajectories
BAND1_VAL1=16.95
BAND1_VAL2=BAND1_VAL1
#BAND1_VAL2=0.05

# Optional band2 
BAND2=True
 # Set True to use band2
BAND2_NAME='longitude' # Dimension to average over: 'latitude' or 'longitude'
BAND2_VAL1='121.75'
BAND2_VAL2=BAND2_VAL1

if SUBDIR=='std':
    FILEPRE='' # e.g., '', '_rac', '_rac_b20_200_n241', '_rac_rm5_n5'
    YEAR_BEG=2001; YEAR_END=2020
    #MONTH1=MONTH2=-999 # Set both MONTH1 and MONTH2 to same (irrelevant) value if outfile_frequency is 'year'
    MONTH1=1; MONTH2=12 # Set month ranges if outfile_frequency is less than 'year'
elif SUBDIR=='processed':
    FILEPRE='_'+TDOMAINID+'_lag' # e.g., TDOMAINID of time mean data.
else:
    raise UserWarning('SUBDIR is invalid.')

CREATE=True

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
descriptor['source']=SOURCE
descriptor['filepre']=FILEPRE
descriptor['band1_name']=BAND1_NAME
descriptor['band1_val1']=BAND1_VAL1
descriptor['band1_val2']=BAND1_VAL2
descriptor['band2']=BAND2
if BAND2:
    descriptor['band2_name']=BAND2_NAME
    descriptor['band2_val1']=BAND2_VAL1
    descriptor['band2_val2']=BAND2_VAL2
    
# Create instance of SpatialSubset object
aa=da.SpatialSubset(**descriptor)

if CREATE:
    if SUBDIR=='std':
        for year in range(YEAR_BEG,YEAR_END+1):
            for month in range(MONTH1,MONTH2+1):
                print('### year={0!s} month={1!s}'.format(year,month))
                aa.year=year
                aa.month=month
                aa.f_spatial_subset()
    elif SUBDIR=='processed':
        aa.f_spatial_subset()
else:
    aa.f_read_spatial_subset()
    aa.data_subset_current=aa.data_subset[-1]

if PLOT:
    x1=aa.data_subset_current
    ndim=x1.data.ndim
    if ndim==2:
        qplt.contourf(x1)
    elif ndim==1:
        qplt.plot(x1)
    
    plt.show()
