"""Get (ERA5) data from Copernicus. 

Run this script from ada.  Can run interactively, but best to run as batch job.

From 20 Aug 2024, moved to use the new CDS-beta Copernicus site to
download ERA5 data, as the old site would no longer be suported from
Sep 2024. Unfortunately the new site has very different metadata in
the netcdf files, so iris will not be able to combine data from the
old and new Copernicus sites. So, created new sources for the ERA5
data from the new site with era5beta replacing just era5 in the
data_source part, e.g., era5beta_trp2_plev_h.

"""

import os

import cdsapi
import cftime
import datetime
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import pdb

import info

BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')
#BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

SDOMAIN='trp2'

VAR_NAME='uwnd'; LEVEL=1000; SOURCE='era5beta'+SDOMAIN+'_plev_h'
#VAR_NAME='vwnd'; LEVEL=1; SOURCE='era5beta'+SDOMAIN+'_sfc_h'
#VAR_NAME='ppt'; LEVEL=1; SOURCE='era5beta'+SDOMAIN+'_sfc_h'

#YEAR_BEG=2017; YEAR_END=2018
YEAR_END=YEAR_BEG+3
MONTH1=1; MONTH2=12 # if outfile_frequency is less than 'year'

DOWNLOAD=True

PLOT=False

#==========================================================================

# Set Copernicus dataset name
if SOURCE in ['era5betatrp_plev_h','era5betaplp_plev_h','era5betamcw_plev_h','era5betaewa_plev_h','era5betaglo_plev_h','era5betauks_plev_h','era5betatrp2_plev_h']:
    dataset='reanalysis-era5-pressure-levels'
elif SOURCE in ['era5betaplp_sfc_h','era5betabar_sfc_h','era5betamcw_sfc_h','era5betaglo_sfc_h']:
    dataset='reanalysis-era5-single-levels'
else:
    raise UserWarning('SOURCE not recognised.')

# Set Copernicus variable name
# To find the Copernicus variable name go to
# cds.climate.copernicus.eu and click on the download tab.
# Select the variable you want then click Show API request.
# This generates code to use below including the Copernicus variable name
if VAR_NAME=='uwnd':
    if SOURCE in ['era5betatrp_plev_h','era5betaplp_plev_h','era5betaewa_plev_h','era5betamcw_plev_h','era5betaglo_plev_h','era5betauks_plev_h','era5betatrp2_plev_h']:
        variable='u_component_of_wind'
    elif SOURCE in ['era5betaplp_sfc_h','era5betabar_sfc_h','era5betamcw_sfc_h','era5betaglo_sfc_h']:
        variable='10m_u_component_of_wind'
    else:
        raise('SOURCE not recognised.')
elif VAR_NAME=='vwnd':
    if SOURCE in ['era5betatrp_plev_h','era5betaplp_plev_h','era5betaewa_plev_h','era5betamcw_plev_h','era5betaglo_plev_h','era5betauks_plev_h','era5betatrp2_plev_h']:
        variable='v_component_of_wind'
    elif SOURCE in ['era5betaplp_sfc_h','era5betabar_sfc_h','era5betamcw_sfc_h','era5betaglo_sfc_h']:
        variable='10m_v_component_of_wind'
    else:
        raise('SOURCE not recognised.')
elif VAR_NAME=='div':
    variable='divergence'
elif VAR_NAME=='vrt':
    variable='vorticity'
elif VAR_NAME=='omega':
    variable='vertical_velocity'
elif VAR_NAME=='ta':
    if SOURCE in ['era5betaplp_sfc_h','era5betabar_sfc_h','era5betamcw_sfc_h']:
        variable='2m_temperature'
    else:
        raise('SOURCE not recognised.')
elif VAR_NAME=='ppt':
    if SOURCE in ['era5betaplp_sfc_h','era5betabar_sfc_h','era5betamcw_sfc_h']:
        variable='mean_total_precipitation_rate'
    else:
        raise('SOURCE not recognised.')
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
        request={
            'product_type': 'reanalysis',
            'data_format': 'netcdf',
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
            request['pressure_level']=[str(LEVEL)]
        print('request: {0!s}'.format(request))
        
        # Set download file name
        target=os.path.join(BASEDIR,SOURCE,'raw',VAR_NAME+'_'+str(LEVEL)+'_'+str(year)+str(month).zfill(2)+'.nc')
        print('target: {0!s}'.format(target))

        if DOWNLOAD:
            # Retrieve data from Copernicus
            c=cdsapi.Client()
            c.retrieve(dataset,request,target)

if PLOT:
    print('# Plot')
    fig=plt.figure()
    x1=iris.load(target)
    x1=x1.concatenate_cube()
    tcoord=x1.coord('time')
    time1=tcoord.units.num2date(tcoord.points[0])
    tconstraint=iris.Constraint(time=time1)
    x2=x1.extract(tconstraint)

    qplt.contourf(x2)
    plt.gca().coastlines()
    plt.show()

    fig.savefig('/gpfs/home/e058/tmp/fig1.png')
