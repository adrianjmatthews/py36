"""Simple plot script. Don't change!"""

import os

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import pdb

BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')
#BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

SOURCE='erainterim_plev_6h'
SUBDIR='std'
VAR_NAME='uwnd'
LEVEL=850


file1=os.path.join(BASEDIR,SOURCE,SUBDIR,VAR_NAME+'_'+str(LEVEL)+'_1998.nc')
print('file1: {0!s}'.format(file1))

cubelist=iris.load(file1)
cube=cubelist.concatenate_cube()

tcoord=cube.coord('time')[0]
time1=tcoord.cell(0)[0]
timecon=iris.Constraint(time=time1)
x2=cube.extract(timecon)
#pdb.set_trace()

fig=plt.figure()

qplt.contourf(x2)
plt.gca().coastlines()
plt.show()

IMAGEFILE=os.path.join(os.path.sep,'gpfs','home','e058','tmp','fig1.png')
fig.savefig(IMAGEFILE)
