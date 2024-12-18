"""Calculate time-lagged mean using data_analysis.TimeDomStats."""

import os

import cftime
import datetime
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import pdb

import data_analysis as da

BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')
#BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

ARCHIVE=True
BASEDIR_ARCHIVE=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

#VAR_NAME='uwnd'; LEVEL=850; SOURCE='era5gloeraiER1_plev_3h'; TDOMAINID='CCER75Elat-15-15-00-22-0.05-00UTC'
#VAR_NAME='uwnd'; LEVEL=850; SOURCE='erainterimEK1_plev_6h'; TDOMAINID='CCEK75Elat-10-10-98-18-0.3-6h'
#VAR_NAME='ppt'; LEVEL=1; SOURCE='imergv07btrmrgp1_sfc_d'; TDOMAINID='CCER75Elat-15-15-00-23-0.05-00UTC'
#VAR_NAME='zg'; LEVEL=200; SOURCE='ncepdoe_plev_d'; TDOMAINID='rmm001a-n2a5'
#VAR_NAME='vwndptap'; LEVEL=850; SOURCE='ncepncar_plev_d'; TDOMAINID='rmm001djf3'
#VAR_NAME='olr'; LEVEL=0; SOURCE='olrinterp_toa_d'; TDOMAINID='rmm001djf3'
#VAR_NAME='ppt'; LEVEL=1; SOURCE='trmm3b42v7p1_sfc_d'; TDOMAINID='CCER'+str(LONC)+'Elat-15-15-98-20-0.05-00UTC'
VAR_NAME='ppt'; LEVEL=1; SOURCE='igcm0002_sfc_3h'; TDOMAINID='CCEK75Elat-5-5-05-19-5.0-3h-360day'

FILEPRE='_rac' # e.g., '', '_rac', '_l30_n241'

# Usually will calculate time mean from e.g., daily data. Set DATA_FROM_ANNCYCLE to False for this.
# To calculate time mean from annual cycle (i.e., to get a mean background state for selected dates)
#   set DATA_FROM_ANNCYCLE to a 2-tuple of (year_beg,year_end) from which annual cycle was calculated
#   e.g., (1998,2018)
DATA_FROM_ANNCYCLE=False
#DATA_FROM_ANNCYCLE=(1998,2018)

METHOD=2
# If using method 1 later (not advised now)
# LAGS is list of integer multiples of timedelta attribute, e.g., 6 hours for erainterim_plev_6h
#LAGS=[0,10,20]
#
# If using method 2, LAGS is a 2- (or 3-) tuple of start timedelta, end timedelta, with
#    optional step timedelta
LAGS=( datetime.timedelta(days=-2), datetime.timedelta(days=2), datetime.timedelta(days=1) )

# Set lazy load to false to force hard load of entire data set.
# Memory hungry but can lead to significant speed up. NB should not need to use this now with method 2
# Need to manually set first and last times to encompass all dates from time domain with lags
LAZY_LOAD=True
if not LAZY_LOAD:
    TIME_FIRST=cftime.DatetimeGregorian(1998,1,1)
    TIME_LAST=cftime.DatetimeGregorian(2020,12,31)

VERBOSE=2

PLOT=False

#==========================================================================

descriptor={}
descriptor['verbose']=VERBOSE
descriptor['var_name']=VAR_NAME
descriptor['level']=LEVEL
descriptor['source']=SOURCE
descriptor['tdomainid']=TDOMAINID
descriptor['basedir']=BASEDIR
descriptor['archive']=ARCHIVE
descriptor['basedir_archive']=BASEDIR_ARCHIVE
descriptor['filepre']=FILEPRE
descriptor['data_from_anncycle']=DATA_FROM_ANNCYCLE
if 'TIME_FIRST' in locals():
    descriptor['time_first']=TIME_FIRST
if 'TIME_LAST' in locals():
    descriptor['time_last']=TIME_LAST

# Create instance of TimeDomStats object
aa=da.TimeDomStats(lazy_load=LAZY_LOAD,**descriptor)

# Calculate event means and time mean
aa.f_lagged_mean(method=METHOD,lags=LAGS)

if PLOT:
    print('# Plot')
    tcoord=aa.lagged_mean.coord('time')
    t1=tcoord.units.num2date(tcoord.points[-1])
    timecon=iris.Constraint(time=t1)
    with iris.FUTURE.context(cell_datetime_objects=True):
        x1=aa.lagged_mean.extract(timecon)

    qplt.contourf(x1); plt.gca().coastlines(); plt.show()
