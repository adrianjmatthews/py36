"""Time filter data using data_analysis.TimeFilter."""

import datetime
import os

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt

import data_analysis as da

BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')
#BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

ARCHIVE=True
BASEDIR_ARCHIVE=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

VAR_NAME='vrt'; LEVEL=850; SOURCE='erainterim_plev_d'
#VAR_NAME='vrt'; LEVEL=850; SOURCE='ncepdoe_plev_d'
#VAR_NAME='ta'; LEVEL=50; SOURCE='ncepncar_plev_d'
#VAR_NAME='wndspd'; LEVEL=1; SOURCE='ncepncar_sfc_d'
#VAR_NAME='olr'; LEVEL=0; SOURCE='olrinterp_toa_d'
#VAR_NAME='sst'; LEVEL=1; SOURCE='sstrey_sfc_d'
#VAR_NAME='ppt'; LEVEL=1; SOURCE='trmm3b42v7_sfc_d'

FILTER='b20_200_n241' # 'rm5_n5' 'h20_n241' 'b20_200_n241' etc.

FILEPRE='_rac' # e.g., '', '_rac',
SUBTRACT=False

YEAR_BEG=2010; YEAR_END=2010
MONTH1=MONTH2=-999 # Set both MONTH1 and MONTH2 to same (irrelevant) value if outfile_frequency is 'year'
#MONTH1=1; MONTH2=1 # Set month ranges if outfile_frequency is less than 'year'

SPLITBLOCK=True

VERBOSE=2

PLOT=False

#==========================================================================

descriptor={}
descriptor['verbose']=VERBOSE
descriptor['file_weights']=os.path.join(os.path.sep,'gpfs','home','e058','home','data','weights','w_'+FILTER+'.txt')
descriptor['var_name']=VAR_NAME
descriptor['level']=LEVEL
descriptor['source']=SOURCE
descriptor['basedir']=BASEDIR
descriptor['archive']=ARCHIVE
descriptor['basedir_archive']=BASEDIR_ARCHIVE
descriptor['filepre']=FILEPRE
descriptor['filter']=FILTER
descriptor['splitblock']=FILTER

# Create instance of TimeFilter object
aa=da.TimeFilter(**descriptor)

# Overwrite irrelevant MONTH1,MONTH2 if outfile_frequency is 'year'
if aa.outfile_frequency=='year':
    MONTH1=MONTH2=-999

for year in range(YEAR_BEG,YEAR_END+1):
    for month in range(MONTH1,MONTH2+1):
        print('### year={0!s} month={1!s}'.format(year,month))
        aa.year=year
        aa.month=month
        aa.time_filter(subtract=SUBTRACT)

if PLOT:
    print('# Plot')
    time_constraint=iris.Constraint(time = lambda cell: aa.timeout1 <= cell <= aa.timeout2)
    tol=0.1
    lon0=0.0
    lon_constraint=iris.Constraint(longitude = lambda cell: lon0-tol <= cell <= lon0+tol)
    lat0=55.0
    lat_constraint=iris.Constraint(latitude = lambda cell: lat0-tol <= cell <= lat0+tol)
    x1=aa.data_in.extract(time_constraint & lon_constraint & lat_constraint)
    x1=x1.concatenate_cube()
    x2=aa.data_out.extract(lon_constraint & lat_constraint)
    qplt.plot(x1,label='in')
    qplt.plot(x2,label='out')
    plt.legend()
    plt.axis('tight')
    qplt.show()
    
