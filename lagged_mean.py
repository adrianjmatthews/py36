"""Calculate time-lagged mean using data_analysis.TimeDomStats."""

import os

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt

import data_analysis as da

BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')
#BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

ARCHIVE=True
BASEDIR_ARCHIVE=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

#VAR_NAME='m_uwndprime_dvrtdxbar'; 
LEVEL=850; SOURCE='erainterimEK2_plev_6h'; 
#TDOMAINID='CCEK75E98-18-0.5-6h'
#VAR_NAME='vrt'; LEVEL=850; SOURCE='ncepdoe_plev_d'; TDOMAINID='cckw75Epm10gt109812vrtsng'
#VAR_NAME='vwndptap'; LEVEL=850; SOURCE='ncepncar_plev_d'; TDOMAINID='rmm001djf3'
#VAR_NAME='olr'; LEVEL=0; SOURCE='olrinterp_toa_d'; TDOMAINID='rmm001djf3'
#VAR_NAME='ppt'; LEVEL=1; SOURCE='trmm3b42v7p1_sfc_d'; TDOMAINID='CCEK75E98-18a-0.5-00UTC'

FILEPRE='' # e.g., '', '_rac', '_l30_n241'

# Usually will calculate time mean from e.g., daily data. Set DATA_FROM_ANNCYCLE to False for this.
# To calculate time mean from annual cycle (i.e., to get a mean background state for selected dates)
#   set DATA_FROM_ANNCYCLE to a 2-tuple of (year_beg,year_end) from which annual cycle was calculated
#   e.g., (1998,2018)
DATA_FROM_ANNCYCLE=False
#DATA_FROM_ANNCYCLE=(1998,2018)

# lags are integer multiples of timedelta attribute, e.g., 6 hours for erainterim_plev_6h
LAGS=[-4,0,4]

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

# Create instance of TimeDomStats object
aa=da.TimeDomStats(**descriptor)

# Calculate event means and time mean
aa.f_lagged_mean(lags=LAGS)

if PLOT:
    print('# Plot')
    tcoord=aa.lagged_mean.coord('time')
    t1=tcoord.units.num2date(tcoord.points[-1])
    timecon=iris.Constraint(time=t1)
    with iris.FUTURE.context(cell_datetime_objects=True):
        x1=aa.lagged_mean.extract(timecon)

    qplt.contourf(x1); plt.gca().coastlines(); plt.show()
