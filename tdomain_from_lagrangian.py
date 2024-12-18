"""Create time domain from CCEW Lagrangian data base.

Optionally, if dbsubset is True, also create subset of CCEW Lagrangian
data base with only the subset of trajectories that are used in the
time domain, and added keywords from TDOMAIN_PARAMS.

"""

import os

import cftime
import datetime
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt

import data_analysis as da

BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')
#BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

VAR_NAME='ppt'; LEVEL=1; SOURCE='igcm0002_sfc_3h'

FILEPRE='_rac' # e.g., '', '_rac', '_rac_minus_l30_n241'

#TIME1=cftime.DatetimeGregorian(1998,1,1)
#TIME2=TIME1+datetime.timedelta(23*365+6-1)-datetime.timedelta(seconds=1) # 30 Dec 2020 TRMM precip 
#TIME2=TIME1+datetime.timedelta(23*365+280-1)-datetime.timedelta(seconds=1) # 30 Sep 2021 TRMM precip 
#
#TIME1=cftime.DatetimeGregorian(2000,7,1)
#TIME2=TIME1+datetime.timedelta(23*365+68-1)-datetime.timedelta(seconds=1) # 31 Aug 2023 imergv07b
#TIME2=TIME1+datetime.timedelta(23*365+190-1)-datetime.timedelta(seconds=1) # 31 Dec 2023 imergv07b
TIME1=cftime.Datetime360Day(3005,1,1) # IGCM
TIME2=cftime.Datetime360Day(3019,12,30)

#LAT1=-2.625; LAT2=-LAT1 # EK
#LAT1=-10; LAT2=-LAT1 # EK
LAT1=-5; LAT2=-LAT1 # EK igcm
#LAT1=-15; LAT2=-LAT1 # ER
#LAT1=-15; LAT2=15 # ER
# Option to use symmetric or antisymmetric component in latitude
BAND1_SYM=False # False, 'sym', or 'antisym'

WAVE_TYPE='EK' # 'EK' or 'ER'

LONC=75
SEASON=False # False (for all year round) or, e.g., 'djf', mam', etc.

#TDOMAIN_PARAMS={'lonc':LONC, 'threshold':0.5, 'threshold_units':'mm hr-1', 'min_lon_extent':False, 'round_to_nearest_time':False, 'dbsubset':True, 'season':SEASON} # EK
#TDOMAIN_PARAMS={'lonc':LONC, 'threshold':False, 'threshold_units':False, 'min_lon_extent':False, 'round_to_nearest_time':'d', 'dbsubset':False}
#TDOMAIN_PARAMS={'lonc':LONC, 'threshold':False, 'threshold_units':False, 'min_lon_extent':False, 'round_to_nearest_time':False, 'dbsubset':False}
#TDOMAIN_PARAMS={'lonc':LONC, 'threshold':0.3/2, 'threshold_units':'mm hr-1', 'min_lon_extent':False, 'round_to_nearest_time':False, 'dbsubset':True, 'season':SEASON} # EK
#TDOMAIN_PARAMS={'lonc':LONC, 'threshold':0.05, 'threshold_units':'mm hr-1', 'min_lon_extent':False, 'round_to_nearest_time':'d', 'dbsubset':True, 'season':SEASON} # ER
TDOMAIN_PARAMS={'lonc':LONC, 'threshold':5.0, 'threshold_units':'mm d-1', 'min_lon_extent':False, 'round_to_nearest_time':False, 'dbsubset':True, 'season':SEASON} # EK igcm

VERBOSE=2

#==========================================================================

descriptor={}
descriptor['verbose']=VERBOSE
descriptor['var_name']=VAR_NAME
descriptor['level']=LEVEL
descriptor['source']=SOURCE
descriptor['basedir']=BASEDIR
descriptor['filepre']=FILEPRE
descriptor['time1']=TIME1
descriptor['time2']=TIME2
descriptor['lat1']=LAT1
descriptor['lat2']=LAT2
descriptor['band1_sym']=BAND1_SYM
descriptor['wave_type']=WAVE_TYPE

# Create instance of CCEWLagrangian object
aa=da.CCEWLagrangian(**descriptor)

# Read Lagrangian data base
aa.f_read_lagrangian()

# Calculate time domain from CCEW events
aa.f_create_time_domain(TDOMAIN_PARAMS)
