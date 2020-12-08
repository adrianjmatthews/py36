import datetime
import os

import cdsapi
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import pdb

import info

BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

SDOMAIN='plp'

#VAR_NAME='vwnd'; LEVEL=925; SOURCE='era5'+SDOMAIN+'_plev_h'
VAR_NAME='vwnd'; LEVEL=1; SOURCE='era5'+SDOMAIN+'_sfc_h'

YEAR=2018
MONTH=9

DOWNLOAD=True

PLOT=True

#==========================================================================

# Set Copernicus dataset name
if SOURCE in ['era5trp_plev_h','era5plp_plev_h']:
    dataset='reanalysis-era5-pressure-levels'
elif SOURCE in ['era5plp_sfc_h']:
    dataset='reanalysis-era5-single-levels'
else:
    raise UserWarning('SOURCE not recognised.')

# Set Copernicus variable name
if VAR_NAME=='uwnd':
    if SOURCE in ['era5trp_plev_h','era5plp_plev_h']:
        variable='u_component_of_wind'
    elif SOURCE in ['era5plp_sfc_h']:
        variable='10m_u_component_of_wind'
elif VAR_NAME=='vwnd':
    if SOURCE in ['era5trp_plev_h','era5plp_plev_h']:
        variable='v_component_of_wind'
    elif SOURCE in ['era5plp_sfc_h']:
        variable='10m_v_component_of_wind'
else:
    raise UserWarning('VAR_NAME not recognised.')

# Create list of day numbers
t1=datetime.datetime(YEAR,MONTH,1)
if MONTH<12:
    t2=datetime.datetime(YEAR,MONTH+1,1)-datetime.timedelta(days=1)
else:
    t2=datetime.datetime(YEAR+1,1,1)-datetime.timedelta(days=1)
ndays=t2.day
print('YEAR,MONTH,ndays: {0!s}, {1!s}, {2!s}'.format(YEAR,MONTH,ndays))
daylist=[str(xx).zfill(2) for xx in range(1,ndays+1)]
print('daylist: {0!s}'.format(daylist))

# Spatial domain for subsetting
lon1=info.sdomains[SDOMAIN]['lon1']
lon2=info.sdomains[SDOMAIN]['lon2']
lat1=info.sdomains[SDOMAIN]['lat1']
lat2=info.sdomains[SDOMAIN]['lat2']

# Create download dictionary
downloaddir={
        'product_type': 'reanalysis',
        'format': 'netcdf',
        'variable': variable,
        'pressure_level': [str(LEVEL)],
        'year': [str(YEAR)],
        'month': str(MONTH).zfill(2),
        'day': daylist,
        'time': [
            '00:00', '01:00', '02:00',
            '03:00', '04:00', '05:00',
            '06:00', '07:00', '08:00',
            '09:00', '10:00', '11:00',
            '12:00', '13:00', '14:00',
            '15:00', '16:00', '17:00',
            '18:00', '19:00', '20:00',
            '21:00', '22:00', '23:00',
        ],
        'area': [lat2, lon1, lat1, lon2],
    }
print('downloaddir: {0!s}'.format(downloaddir))

# Set download file name
filei1=os.path.join(BASEDIR,SOURCE,'raw',VAR_NAME+'_'+str(LEVEL)+'_'+str(YEAR)+str(MONTH).zfill(2)+'.nc')
print('filei1: {0!s}'.format(filei1))

if DOWNLOAD:
    # Retrieve data from Copernicus
    c=cdsapi.Client()
    c.retrieve(dataset,downloaddir,filei1)

if PLOT:
    print('# Plot')
    fig=plt.figure()
    x1=iris.load(filei1)
    x1=x1.concatenate_cube()
    tcoord=x1.coord('time')
    time1=tcoord.units.num2date(tcoord.points[0])
    tconstraint=iris.Constraint(time=time1)
    x2=x1.extract(tconstraint)

    qplt.contourf(x2)
    plt.gca().coastlines()
    plt.show()

    fig.savefig('/gpfs/home/e058/tmp/fig1.png')
