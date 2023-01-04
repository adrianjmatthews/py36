"""Get (ERA5) data from Copernicus. 

Run this script from ada.  Can run interactively, but best to run as batch job.

"""

import os

import cdsapi
import cftime
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import pdb

import info

BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')
#BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

SDOMAIN='ewa'

VAR_NAME='vwnd'; LEVEL=850; SOURCE='era5'+SDOMAIN+'_plev_h'
#VAR_NAME='ta'; LEVEL=1; SOURCE='era5'+SDOMAIN+'_sfc_h'
#VAR_NAME='ppt'; LEVEL=1; SOURCE='era5'+SDOMAIN+'_sfc_h'

YEAR_BEG=2020; YEAR_END=2020 # if outfile_frequency is 'year' or less
MONTH1=2; MONTH2=5 # if outfile_frequency is less than 'year'

DOWNLOAD=True

PLOT=False

#==========================================================================

# Set Copernicus dataset name
if SOURCE in ['era5trp_plev_h','era5plp_plev_h','era5mcw_plev_h','era5ewa_plev_h']:
    dataset='reanalysis-era5-pressure-levels'
elif SOURCE in ['era5plp_sfc_h','era5bar_sfc_h','era5mcw_sfc_h']:
    dataset='reanalysis-era5-single-levels'
else:
    raise UserWarning('SOURCE not recognised.')

# Set Copernicus variable name
# To find the Copernicus variable name go to
# cds.climate.copernicus.eu and click on the download tab.
# Select the variable you want then click Show API request.
# This generates code to use below including the Copernicus variable name
if VAR_NAME=='uwnd':
    if SOURCE in ['era5trp_plev_h','era5plp_plev_h','era5ewa_plev_h']:
        variable='u_component_of_wind'
    elif SOURCE in ['era5plp_sfc_h','era5bar_sfc_h','era5mcw_sfc_h']:
        variable='10m_u_component_of_wind'
elif VAR_NAME=='vwnd':
    if SOURCE in ['era5trp_plev_h','era5plp_plev_h','era5ewa_plev_h']:
        variable='v_component_of_wind'
    elif SOURCE in ['era5plp_sfc_h','era5bar_sfc_h','era5mcw_sfc_h']:
        variable='10m_v_component_of_wind'
elif VAR_NAME=='div':
    variable='divergence'
elif VAR_NAME=='ta':
    if SOURCE in ['era5plp_sfc_h','era5bar_sfc_h','era5mcw_sfc_h']:
        variable='2m_temperature'
elif VAR_NAME=='ppt':
    if SOURCE in ['era5plp_sfc_h','era5bar_sfc_h','era5mcw_sfc_h']:
        variable='mean_total_precipitation_rate'
else:
    raise UserWarning('VAR_NAME not recognised.')

# Spatial domain for subsetting
lon1=info.sdomains[SDOMAIN]['lon1']
lon2=info.sdomains[SDOMAIN]['lon2']
lat1=info.sdomains[SDOMAIN]['lat1']
lat2=info.sdomains[SDOMAIN]['lat2']

# Loop over years and months
for year in range(YEAR_BEG,YEAR_END+1):
    for month in range(MONTH1,MONTH2+1):
        print('### year={0!s} month={1!s}'.format(year,month))

        # Create list of day numbers
        t1=cftime.DatetimeGregorian(year,month,1)
        if month<12:
            t2=cftime.DatetimeGregorian(year,month+1,1)-datetime.timedelta(days=1)
        else:
            t2=cftime.DatetimeGregorian(year+1,1,1)-datetime.timedelta(days=1)
        ndays=t2.day
        print('year,month,ndays: {0!s}, {1!s}, {2!s}'.format(year,month,ndays))
        daylist=[str(xx).zfill(2) for xx in range(1,ndays+1)]
        print('daylist: {0!s}'.format(daylist))
        
        # Create download dictionary
        downloaddir={
            'product_type': 'reanalysis',
            'format': 'netcdf',
            'variable': variable,
            'year': [str(year)],
            'month': str(month).zfill(2),
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
        if SOURCE.split('_')[1]=='plev':
            downloaddir['pressure_level']=[str(LEVEL)]
        print('downloaddir: {0!s}'.format(downloaddir))
        
        # Set download file name
        filei1=os.path.join(BASEDIR,SOURCE,'raw',VAR_NAME+'_'+str(LEVEL)+'_'+str(year)+str(month).zfill(2)+'.nc')
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
