"""Create time domain from CCEW Lagrangian data base."""

import datetime
import os

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt

import data_analysis as da

BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

VAR_NAME='ppt'; LEVEL=1; SOURCE='trmm3b42v7p2_sfc_3h'

FILEPRE='' # e.g., '', '_rac', '_rac_minus_l30_n241'

TIME1=datetime.datetime(1998,1,1)
TIME2=TIME1+datetime.timedelta(23*365+6-1)-datetime.timedelta(seconds=1)

LAT1=-2.625; LAT2=-LAT1 # EK
LAT1=-15; LAT2=-LAT1 # ER

WAVE_TYPE='ER' # 'EK' or 'ER'

#TDOMAIN_PARAMS={'lonc':75, 'threshold':0.5, 'threshold_units':'mm hr-1', 'min_lon_extent':False, 'round_to_nearest_time':'6h'}
#TDOMAIN_PARAMS={'lonc':102, 'threshold':False, 'threshold_units':False, 'min_lon_extent':False, 'round_to_nearest_time':'d'}
#TDOMAIN_PARAMS={'lonc':106, 'threshold':False, 'threshold_units':False, 'min_lon_extent':False, 'round_to_nearest_time':False}
#TDOMAIN_PARAMS={'lonc':106, 'threshold':0.4, 'threshold_units':'mm hr-1', 'min_lon_extent':False, 'round_to_nearest_time':'d'}
TDOMAIN_PARAMS={'lonc':75, 'threshold':0.05, 'threshold_units':'mm hr-1', 'min_lon_extent':False, 'round_to_nearest_time':'d'} # ER

VERBOSE=2

PLOT=False

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
descriptor['wave_type']=WAVE_TYPE

# Create instance of CCEWLagrangian object
aa=da.CCEWLagrangian(**descriptor)

# Read Lagrangian data base
aa.f_read_lagrangian()

# Calculate time domain from CCEW events
aa.f_create_time_domain(TDOMAIN_PARAMS)
