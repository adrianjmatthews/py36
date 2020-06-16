"""Calculate time mean statistics using data_analysis.TimeDomStats."""

import datetime
import os

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt

import data_analysis as da

BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')
#BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

ARCHIVE=True
BASEDIR_ARCHIVE=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

#VAR_NAME='uwnd'; 
LEVEL=500; SOURCE='erainterim_plev_d'; 
#TDOMAINID='ann9818'
#VAR_NAME='swpd'; LEVEL='all'; SOURCE='sg613m031oi01_zlev_h'; TDOMAINID='boballsg'
#VAR_NAME='ppt'; LEVEL=1; SOURCE='trmm3b42v7p1_sfc_d'; TDOMAINID='jjas98-12'
#VAR_NAME='vwnd'; LEVEL=850; SOURCE='erainterim_plev_d'; TDOMAINID='cckw75Epm10gt109812'
#VAR_NAME='zg'; LEVEL=250; SOURCE='hadgem2esajhog_plev_d'; TDOMAINID='jja86-88-360day'
#VAR_NAME='zg'; LEVEL=250; SOURCE='ncepncar_plev_d'; TDOMAINID='ann9917'
#VAR_NAME='uwnd'; LEVEL=1; SOURCE='ncepncar_sfc_d'; TDOMAINID='apr9815'
#VAR_NAME='olr'; LEVEL=0; SOURCE='olrinterp_toa_d'; TDOMAINID='rmm004djf3'
#VAR_NAME='sst'; LEVEL=1; SOURCE='sstrey_sfc_d'; TDOMAINID='djf8182-1516'

FILEPRE='' # e.g., '', '_rac', '_rac_minus_l30_n241'

# Optional parameters for use in null distribution calculation
#NMC=10; PERCENTILES_NULL=[1,2.5,5,10,20,30,50,70,80,90,95,97.5,99]; MAX_DAY_SHIFT=15
#TIME_FIRST=datetime.datetime(1998,1,1)
#TIME_LAST=datetime.datetime(2012,12,31)

LAZY_LOAD=True
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
if 'NMC' in locals():
    descriptor['nmc']=NMC
if 'PERCENTILES_NULL' in locals():
    descriptor['percentiles_null']=PERCENTILES_NULL
if 'MAX_DAY_SHIFT' in locals():
    descriptor['max_day_shift']=MAX_DAY_SHIFT
if 'TIME_FIRST' in locals():
    descriptor['time_first']=TIME_FIRST
if 'TIME_LAST' in locals():
    descriptor['time_last']=TIME_LAST
if 'NODE_NUMBER' in locals():
    descriptor['node_number']=NODE_NUMBER

# Create instance of TimeDomStats object
aa=da.TimeDomStats(lazy_load=LAZY_LOAD,**descriptor)

# Calculate event means and time mean
aa.event_means()
aa.f_time_mean()

# Calculate components of null distribution (use with batch distribution)
#aa.f_time_mean_null_distribution_component()

# Calculate percentiles of null distribution
#aa.f_percentiles_null()

if PLOT:
    print('# Plot')
    fig=plt.figure()
    x1=aa.time_mean

    #per_constraint=iris.Constraint(percentile=2.5)
    #x1=aa.mean_percentiles_null.extract(per_constraint)

    qplt.contourf(x1)
    plt.gca().coastlines()
    plt.show()

    fig.savefig('/gpfs/home/e058/tmp/fig1.png')
