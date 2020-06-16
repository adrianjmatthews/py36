"""Calculate time-lagged mean using data_analysis.TimeDomStats."""

import os

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt

import data_analysis as da

BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')

ARCHIVE=True
BASEDIR_ARCHIVE=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

#VAR_NAME='uwnd'; 
LEVEL=500; SOURCE='erainterimEK2_plev_6h'; TDOMAINID='CCEK75E98-18-0.5-6h'
#VAR_NAME='vrt'; LEVEL=850; SOURCE='ncepdoe_plev_d'; TDOMAINID='cckw75Epm10gt109812vrtsng'
#VAR_NAME='vwndptap'; LEVEL=850; SOURCE='ncepncar_plev_d'; TDOMAINID='rmm001djf3'
#VAR_NAME='olr'; LEVEL=0; SOURCE='olrinterp_toa_d'; TDOMAINID='rmm001djf3'

FILEPRE='' # e.g., '', '_rac',

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
