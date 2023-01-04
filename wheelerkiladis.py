"""Calculate Wheeler-Kiladis diagnostics using data_analysis.WheelerKiladis."""

import os

import cftime
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import pdb

import data_analysis as da

BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')
#BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

ARCHIVE=True
BASEDIR_ARCHIVE=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

VAR_NAME='ppt'; LEVEL=1; SOURCE='trmm3b42v7p1_sfc_3h'
#VAR_NAME='uwnd'; LEVEL=850; SOURCE='erainterim_plev_6h'

FILEPRE='' # e.g., '', '_rac', '_rac_minus_l30_n241'

TIME1=cftime.DatetimeGregorian(1998,1,1)
#TIME2=TIME1+datetime.timedelta(21*365+6-1+272)-datetime.timedelta(seconds=1) # 29 Sep 2019 TRMM precip 
TIME2=TIME1+datetime.timedelta(23*365+6-1)-datetime.timedelta(seconds=1) # 30 Dec 2020 TRMM precip 
#TIME2=TIME1+datetime.timedelta(21*365+6-1)-datetime.timedelta(seconds=1) # erainterim

#LAT1=-2.625; LAT2=-LAT1
#LAT1=-10; LAT2=-LAT1
#LAT1=-15; LAT2=-LAT1
LAT1=-15; LAT2=15
#LAT1='-00.3508'
#LAT1=2.625; LAT2=7.375
#LAT2=LAT1
# Option to use symmetric or antisymmetric component in latitude
BAND1_SYM='sym' # False, 'sym', or 'antisym'

HOVMOLLER_PARAMS={'cosine_tapering':False, 'cosine_tapering_fraction':0.1}

#WAVE_TYPE='none'; WAVE_PARAMS={}

#WAVE_TYPE='EK'; WAVE_PARAMS={'cphasex_min':5.0, 'cphasex_max':30.0, 'cphasex_units':'m s-1', 'ss_min':1, 'ss_max':14, 'freq_min':0.0333333/8, 'freq_max':0.4/8, 'freq_units':'cycles per 3h'} # TRMM precip
WAVE_TYPE='ER'; WAVE_PARAMS={'nn':1, 'H_min':2.5, 'H_max':90.0, 'H_units':'m', 'ss_min':-10, 'ss_max':-1, 'freq_min':0.01666667/8, 'freq_max':0.4/8, 'freq_units':'cycles per 3h'} # TRMM precip

#WAVE_TYPE='EK'; WAVE_PARAMS={'cphasex_min':5.0, 'cphasex_max':30.0, 'cphasex_units':'m s-1', 'ss_min':1, 'ss_max':14, 'freq_min':0.0333333/4, 'freq_max':0.4/4, 'freq_units':'cycles per 6h'} # erainterim
#WAVE_TYPE='ER'; WAVE_PARAMS={'nn':1, 'H_min':2.5, 'H_max':90.0, 'H_units':'m', 'ss_min':-10, 'ss_max':-1, 'freq_min':0.01666667/4, 'freq_max':0.4/4, 'freq_units':'cycles per 6h'} # erainterim

PLOT=False

VERBOSE=2

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
descriptor['band1_sym']=BAND1_SYM
descriptor['hovmoller_params']=HOVMOLLER_PARAMS
descriptor['wave_type']=WAVE_TYPE
descriptor['wave_params']=WAVE_PARAMS

# Create instance of WheelerKiladis object
aa=da.WheelerKiladis(**descriptor)

# 'calculate' and write 2-D Hovmoller as single cube, or 'read' pre-existing,
# or create 'analytical' test Hovmoller
aa.f_hovmoller(mode='calculate')

# Calculate 2-D FFT
aa.f_fft()

# Filter for selected wave
aa.f_filter()

# Calculate inverse 2-D FFT of wave-filtered data
aa.f_ifft()

if PLOT:
    print('# Plot')

    #x1=aa.data_hov
    #x1=aa.data_hovfftWKmask
    #x1=aa.cphasex2d
    #x2=aa.data_hovWKfilt
    #x3=aa.data_hovdiff

    harm_constraint=iris.Constraint(integer_harmonic_in_time=lambda cell: cell<=20)
    freq_constraint=iris.Constraint(frequency=lambda cell: cell<=0.5)
    ss_constraint=iris.Constraint(integer_zonal_wavenumber=lambda cell: -25<=cell<=25)
    time_constraint=iris.Constraint(time=lambda cell: cftime.DatetimeGregorian(2011,11,11)<=cell.point<=cftime.DatetimeGregorian(2011,12,21))
    lon_constraint=iris.Constraint(longitude=lambda cell: 40<=cell<=120)

    x1=aa.data_hov.extract(time_constraint & lon_constraint)
    x2=aa.data_hovWKfilt.extract(time_constraint & lon_constraint)
    #x3=aa.data_hovdiff.extract(time_constraint & lon_constraint)

    kount=1
    for xx in [x1,x2]:
        print('kount: {0!s}'.format(kount))
        fig=plt.figure()
        qplt.contourf(xx)
        plt.show()
        fig.savefig('/gpfs/home/e058/tmp/fig'+str(kount)+'.png')
        kount+=1

print('Successfully completed.')
