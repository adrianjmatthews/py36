"""Average data to lower longitude resolution using .

This is a bit of hack code just to regrid longitude in a single cube
Hovmoller.

"""

import os

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt

import data_analysis as da

BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')
#BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

VAR_NAME='ppt'; LEVEL=1; SOURCE1='trmm3b42v7p1_sfc_3h'; SOURCE2='trmm3b42v7p2_sfc_3h'

#FILEPRE='hovWKfiltEK_lat_-2.625_2.625_1998-01-01_2021-09-30'
FILEPRE='hov_lat_-2.625_2.625_1998-01-01_2021-09-30'

#FILEPRE='hovWKfiltEK_lat_-10_10_1998-01-01_2020-12-30'
#FILEPRE='hov_lat_-10_10_1998-01-01_2020-12-30'

#FILEPRE='hovWKfiltER_lat_-15_15_sym_1998-01-01_2020-12-30'
#FILEPRE='hov_lat_-15_15_sym_1998-01-01_2020-12-30'

NAVE=4

PLOT=False

VERBOSE=2

#------------------------------------------------------------------

filei1=os.path.join(BASEDIR,SOURCE1,'processed',VAR_NAME+'_'+str(LEVEL)+'_'+FILEPRE+'.nc')
print('filei1: {0!s}'.format(filei1))
x1=iris.load_cube(filei1)
x2=da.f_longitude_average(x1,NAVE)
fileo1=os.path.join(BASEDIR,SOURCE2,'processed',VAR_NAME+'_'+str(LEVEL)+'_'+FILEPRE+'.nc')
print('fileo1: {0!s}'.format(fileo1))
iris.save(x2,fileo1)

if PLOT:
    print('# Plot')
    fig=plt.figure()

    qplt.contourf(x1)
    plt.show()

    fig.savefig('/gpfs/home/e058/tmp/fig1.png')
