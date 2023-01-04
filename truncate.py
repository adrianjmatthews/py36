"""Spectrally truncate iris cube."""

import os

import cftime
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import pdb
from windspharm.iris import VectorWind

import data_analysis as da

BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

TRUNCATION=42 # Triangular truncation limit

SOURCE1='erainterimNEK1_plev_6h'
SUBDIR='processed'
VAR_NAME='vwnd'
TDOMAINID='CCEK75E98-18-0.5-6h'
FILENAME=VAR_NAME+'_850_'+TDOMAINID+'_lag.nc'

PLOT=False

VERBOSE=2

#------------------------------------------------------------------

print('# Calculate source2 from SOURCE1 and TRUNCATION')
xx=SOURCE1.split('_')
source2=xx[0]+'T'+str(TRUNCATION)+'_'+xx[1]+'_'+xx[2]
print('SOURCE1: {0!s}'.format(SOURCE1))
print('TRUNCATION: {0!s}'.format(TRUNCATION))
print('source2: {0!s}'.format(source2))

print('# Read input cube')
filein1=os.path.join(BASEDIR,SOURCE1,SUBDIR,FILENAME)
print('filein1: {0!s}'.format(filein1))
x1=iris.load(filein1)
cube_in=x1[0]

print('# Truncate input cube')
cube_out=da.truncate(cube_in,TRUNCATION)

print('# Save output cube')
fileout1=os.path.join(BASEDIR,source2,SUBDIR,FILENAME)
print('fileout1: {0!s}'.format(fileout1))
iris.save(cube_out,fileout1)

if PLOT:
    print('# Plot')

    #time1=cftime.DatetimeGregorian(1000,1,1) # lag 0 is 1 Jan 1000
    #time_constraint=iris.Constraint(time=time1)
    #x1=cube_in.intersection(longitude=(75,255),latitude=(-50,5))
    #y1=cube_out.intersection(longitude=(75,255),latitude=(-50,5))
    #with iris.FUTURE.context(cell_datetime_objects=True):
    #    x2=x1.extract(time_constraint)
    #    y2=y1.extract(time_constraint)
     
    
    #qplt.contourf(x2)
    #qplt.contour(y2)
    #plt.gca().coastlines()

    qplt.contourf(cube_out)
    plt.gca().coastlines()
    
    plt.show()
    
