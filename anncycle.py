"""Calculate and subtract annual cycle using data_analysis.AnnualCycle."""

import datetime
import os

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt

import data_analysis as da

#BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')
BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')

ARCHIVE=False
BASEDIR_ARCHIVE=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

# YEAR1-YEAR2 are complete years over which to calculate annual cycle
VAR_NAME='vrt'; LEVEL=500; SOURCE='erainterim_plev_6h'; YEAR1=1998; YEAR2=2018
#VAR_NAME='vwnd'; LEVEL=850; SOURCE='erainterim_plev_d'; YEAR1=1998; YEAR2=2018
#VAR_NAME='zg'; LEVEL=250; SOURCE='hadgem2esajhog_plev_d'; YEAR1=1985; YEAR2=1993
#VAR_NAME='vrt'; LEVEL=850; SOURCE='ncepdoe_plev_d'; YEAR1=1980; YEAR2=2015
#VAR_NAME='wndspd'; LEVEL=1; SOURCE='ncepncar_sfc_d'; YEAR1=1979; YEAR2=2015
#VAR_NAME='zg'; LEVEL=250; SOURCE='ncepncar_plev_d'; YEAR1=1979; YEAR2=2017
#VAR_NAME='olr'; LEVEL=0; SOURCE='olrinterp_toa_d'; YEAR1=2010; YEAR2=2011
#VAR_NAME='sst'; LEVEL=1; SOURCE='sstrey_sfc_d'; YEAR1=1982; YEAR2=2015
#VAR_NAME='ppt'; LEVEL=1; SOURCE='trmm3b42v7p1_sfc_3h'; YEAR1=1998; YEAR2=1999

# If TIME1 and TIME2 are False, anomalies from annual cycle are calculated
# over whole input data set,
# If TIME1 is datetime object and TIME2 is False, calculated from TIME1 to end
# If TIME1 is False and TIME2 is datetime, calculated from start to TIME2
# If TIME1 and TIME2 are datetime, calculated from TIME1 to TIME2
TIME1=TIME2=False
#TIME1=datetime.datetime(2017,1,1); TIME2=False
#TIME1=False;  TIME2=datetime.datetime(1980,12,31)
#TIME1=datetime.datetime(1998,1,1);  TIME2=datetime.datetime(1998,12,31)

# ANNCYCLE_SOURCE is e.g., 'erainterim_plev_d' if subtracting annual cycle of daily 
# data from higher frequency (e.g., 6h) data.
# Set to False to switch off
ANNCYCLE_SOURCE='erainterim_plev_d'

NHARM=3
if VAR_NAME=='ppt':
    NHARM=6
print('NHARM',NHARM)

KMIN=5
DETREND=False

VERBOSE=2

PLOT=False

#------------------------------------------------------------------

descriptor={}
descriptor['verbose']=VERBOSE
descriptor['var_name']=VAR_NAME
descriptor['level']=LEVEL
descriptor['source']=SOURCE
descriptor['year1']=YEAR1
descriptor['year2']=YEAR2
descriptor['time1']=TIME1
descriptor['time2']=TIME2
descriptor['anncycle_source']=ANNCYCLE_SOURCE
descriptor['nharm']=NHARM
descriptor['basedir']=BASEDIR
descriptor['kmin']=KMIN
descriptor['detrend']=DETREND
descriptor['archive']=ARCHIVE
descriptor['basedir_archive']=BASEDIR_ARCHIVE

# Create instance of AnnualCycle object
aa=da.AnnualCycle(**descriptor)

if aa.detrend:
    # Either calculate linear trend for later removal, or read previously calculated one
    aa.f_trend()
    #aa.f_read_trend()

# Either create and save raw annual cycle or read in previously calculated one
#aa.f_anncycle_raw()
#aa.f_read_anncycle_raw()

# Either create and save smoothed annual cycle or read in previously calculated one
#aa.f_anncycle_smooth()
aa.f_read_anncycle_smooth()

# Create or read higher frequency smoothed annual cycle from daily annual cycle if needed
if aa.anncycle_source:
    aa.f_expand_anncycle_smooth()
    #aa.f_read_expanded_anncycle_smooth()
    pass

# Either create and save anomaly data (smooothed annual cycle subtracted),
# or read in previously calculated anomaly data
#aa.f_subtract_anncycle()
#aa.f_read_subtract_anncycle()

if PLOT:
    #timecon=da.set_time_constraint(datetime.datetime(1,1,16),False,calendar=aa.calendar,verbose=True)
    timecon=da.set_time_constraint(datetime.datetime(2017,1,16),False,calendar=aa.calendar,verbose=True)
    #timecon=da.set_time_constraint(datetime.datetime(1980,1,1),datetime.datetime(2010,12,31))
    latcon=iris.Constraint(latitude=0.)
    loncon=iris.Constraint(longitude=90.)

    #x1=aa.data_anncycle_raw.extract(timecon)
    #
    #x1=aa.data_anncycle_smooth.extract(timecon)
    #
    x1=aa.data_anncycle_rm.extract(timecon)
    x1=x1[0]
    #
    #x1=aa.cc
    #
    qplt.contourf(x1)
    plt.gca().coastlines()

    #x1=aa.data_anncycle_smooth.extract(latcon & loncon)
    #x2=aa.data_anncycle_raw.extract(latcon & loncon)
    #
    #x1=aa.data_anncycle_rm.extract(timecon & latcon & loncon)
    #x1=x1.concatenate_cube()
    #x2=aa.data_in.extract(timecon & latcon & loncon)
    #x2=x2.concatenate_cube()
    #
    #qplt.plot(x1)
    #qplt.plot(x2)
    
    plt.show()
