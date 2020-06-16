"""Create Lagrangian data base of convectively coupled equatorial waves."""

import datetime
import os

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt

import data_analysis as da

#BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')
BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

ARCHIVE=False
BASEDIR_ARCHIVE=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

VAR_NAME='ppt'; LEVEL=1; SOURCE='trmm3b42v7p2_sfc_3h'

FILEPRE='' # e.g., '', '_rac', '_rac_minus_l30_n241'

TIME1=datetime.datetime(1998,1,1)
TIME2=TIME1+datetime.timedelta(21*365+278-1)-datetime.timedelta(seconds=1)

LAT1=-2.625; LAT2=-LAT1

WAVE_TYPE='EK'; 
EVENT_PARAMS={'threshold':0.15 , 'threshold_units':'mm hr-1' , 'traj_min_time_length':12 , 'traj_min_time_length_units':'h' , 'prune_threshold1':0.5 , 'prune_threshold2':0.25}
#EVENT_PARAMS={'threshold':0.25 , 'threshold_units':'mm hr-1' , 'traj_min_time_length':12 , 'traj_min_time_length_units':'h' , 'prune_threshold1':0.7 , 'prune_threshold2':0.35}

VERBOSE=2

PLOT=False

#==========================================================================

descriptor={}
descriptor['verbose']=VERBOSE
descriptor['var_name']=VAR_NAME
descriptor['level']=LEVEL
descriptor['source']=SOURCE
descriptor['basedir']=BASEDIR
descriptor['archive']=ARCHIVE
descriptor['basedir_archive']=BASEDIR_ARCHIVE
descriptor['filepre']=FILEPRE
descriptor['time1']=TIME1
descriptor['time2']=TIME2
descriptor['lat1']=LAT1
descriptor['lat2']=LAT2
descriptor['wave_type']=WAVE_TYPE

# Create instance of CCEWLagrangian object
aa=da.CCEWLagrangian(**descriptor)

# Read WK filtered hovmoller
aa.f_read_hovWKfilt()

# Create Lagrangian database of CCEW events
aa.f_events(EVENT_PARAMS)

if PLOT:
    print('# Plot')
    fig=plt.figure()

    #x1=aa.data_hovWKfilt
    #x1=aa.data_hovmax

    time_constraint=iris.Constraint(time=lambda cell: datetime.datetime(2011,11,11)<=cell<=datetime.datetime(2011,12,21))
    lon_constraint=iris.Constraint(longitude=lambda cell: 40<=cell<=120)

    x1=aa.data_hovmax.extract(time_constraint & lon_constraint)

    qplt.contourf(x1)
    plt.show()

    fig.savefig('/gpfs/home/e058/tmp/fig1.png')
