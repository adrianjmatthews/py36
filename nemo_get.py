"""Get NEMO reanalysis data. 

Run this script from ada.  Can run interactively, but best to run as batch job.

"""

import os

import cftime
import datetime
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import pdb

BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')
#BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

SOURCE='glorys12v1eq1_zlev_d'

VAR_NAME='vcur'; LEVEL=0.494

YEAR_BEG=2017; YEAR_END=2017
MONTH1=1; MONTH2=1 # if outfile_frequency is less than 'year'

DOWNLOAD=True

PLOT=False

#==========================================================================

# Set variable name
if VAR_NAME=='swtheta':
    variable='thetao'
elif VAR_NAME=='swsal':
    variable='so'
elif VAR_NAME=='ucur':
    variable='uo'
elif VAR_NAME=='vcur':
    variable='vo'
else:
    raise UserWarning('VAR_NAME not recognised.')

# Set level
# Check level is to 3 decimal places exactly (my convention!)
xx=str(LEVEL).split('.')
if len(xx)!=2 or len(xx[1])!=3:
    raise UserWarning('LEVEL should be a decimal with 3 decimal places.')
delta_level=0.001
depth_min=LEVEL-delta_level
depth_max=LEVEL+delta_level

# Loop over years and months
for year in range(YEAR_BEG,YEAR_END+1):
    for month in range(MONTH1,MONTH2+1):
        print('### year={0!s} month={1!s}'.format(year,month))
        date_min=cftime.DatetimeGregorian(year,month,1)
        if month<12:
            date_max=cftime.DatetimeGregorian(year,month+1,1)-datetime.timedelta(seconds=1)
        else:
            date_max=cftime.DatetimeGregorian(year+1,1,1)-datetime.timedelta(seconds=1)
        date_min=str(date_min)
        date_max=str(date_max)
        print('date_min,date_max: {0!s}, {1!s}'.format(date_min,date_max,))

        # Set download file name
        out_dir=os.path.join(BASEDIR,SOURCE,'raw')
        out_name=os.path.join(VAR_NAME+'_'+str(LEVEL)+'_'+str(year)+str(month).zfill(2)+'.nc')
        print('out_dir: {0!s}'.format(out_dir))
        print('out_name: {0!s}'.format(out_name))

        if DOWNLOAD:
            # Retrieve data from CMEMS
            command1="python -m motuclient --motu https://my.cmems-du.eu/motu-web/Motu --service-id GLOBAL_MULTIYEAR_PHY_001_030-TDS --product-id cmems_mod_glo_phy_my_0.083_P1D-m --longitude-min -180 --longitude-max 179.9167 --latitude-min -15 --latitude-max 15 --date-min '" + date_min + "' --date-max '" + date_max + "' --depth-min " + str(depth_min) + " --depth-max " + str(depth_max) + " --variable " + variable + " --out-dir '" + out_dir  + "' --out-name '" + out_name + "'  --user amatthews4 --pwd  'CMEMS_MATTHEWS_2021'"
            print(command1)
            #os.system(command1)

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
