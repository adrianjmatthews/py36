"""Calculate vorticity budget terms."""

import os

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt

import data_analysis as da

BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')

ARCHIVE=False
BASEDIR_ARCHIVE=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

SOURCE='igcm0002_plev_3h'
LEVEL_BELOW=875; LEVEL=850; LEVEL_ABOVE=825

FILEPRE='' # '', '_rac', etc.

FILEPREANNCYCLE=False
#FILEPREANNCYCLE='_ac_smooth_expanded_1998_2018' # if decomposing into annual cycle and perturbation

# Only use SOURCE2 if decomposing into filtered and remainder. Use SOURCE2 for the remainder.
SOURCE2=False
#SOURCE2='era5gloeraiNER1_plev_3h'

YEAR=3005
#YEAR=range(3005,3049+1)

MONTH=[-999] # Dummy value if outfile_frequency is 'year'
#MONTH=1
#MONTH=range(1,12+1) # If outfile_frequency is less than 'year' 

PLOT=False

VERBOSE=2

#------------------------------------------------------------------

descriptor={}
descriptor['verbose']=VERBOSE
descriptor['basedir']=BASEDIR
descriptor['basedir_archive']=BASEDIR_ARCHIVE
descriptor['archive']=ARCHIVE
descriptor['source']=SOURCE
descriptor['level']=LEVEL
descriptor['filepre']=FILEPRE
descriptor['filepreanncycle']=FILEPREANNCYCLE
descriptor['source2']=SOURCE2

# Create instance of CubeDiagnostics object
aa=da.CubeDiagnostics(**descriptor)

# Lazy read data (select which is needed)
aa.f_read_data('uwnd',LEVEL_BELOW)
aa.f_read_data('uwnd',LEVEL)
aa.f_read_data('uwnd',LEVEL_ABOVE)
aa.f_read_data('vwnd',LEVEL_BELOW)
aa.f_read_data('vwnd',LEVEL)
aa.f_read_data('vwnd',LEVEL_ABOVE)
aa.f_read_data('vrt',LEVEL_BELOW)
aa.f_read_data('vrt',LEVEL)
aa.f_read_data('vrt',LEVEL_ABOVE)
aa.f_read_data('omega',LEVEL)
aa.f_read_data('div',LEVEL)

# Read annual cycle data if needed
#aa.f_read_anncycle('vwnd',LEVEL,verbose=VERBOSE)
#aa.f_read_anncycle('vrt',LEVEL,verbose=VERBOSE)
#aa.f_read_anncycle('div',LEVEL,verbose=VERBOSE)

# Read remainder data if needed
#aa.f_read_data_source2('vrt',LEVEL,verbose=VERBOSE)
#aa.f_read_data_source2('div',LEVEL,verbose=VERBOSE)
#aa.f_read_data_source2('uwnd',LEVEL,verbose=VERBOSE)
#aa.f_read_data_source2('vwnd',LEVEL,verbose=VERBOSE)

iter_year=da.iter_generator(YEAR)
iter_month=da.iter_generator(MONTH)
for year in iter_year:
    for month in iter_month:
        print('### year={0!s} month={1!s}'.format(year,month))
        aa.year=year
        aa.month=month
        #
        # Calculate all vorticity budget terms (usual function to call)
        aa.f_vrtbudget(LEVEL_BELOW,LEVEL,LEVEL_ABOVE,flag_hypdiff=True,ndel=6,kappa=3e8)
        #
        # Only calculate m_vrt_div term, for diagnostic purposes
        #aa.f_m_vrt_div(LEVEL)
        #
        # Calculate dvwnddx and m_duwnddy contributions to vorticity, for diagnostic purposes
        #aa.f_vrt_components(LEVEL)
        #
        # Calculate duwnddx and dvwnddy contributions to divergence, for diagnostic purposes
        #aa.f_div_components(LEVEL)
        #
        # Calculate individual terms decomposed into annual cycle and perturbation combinations
        #aa.f_m_vrt_div_annpert(LEVEL)
        #aa.f_m_uwnd_dvrtdx_annpert(LEVEL)
        #aa.f_m_vwnd_dvrtdy_annpert(LEVEL)
        #
        # Calculate individual terms decomposed into perturbation and remainder combinations
        #aa.f_m_vrt_div_pertremainder(LEVEL)
        #aa.f_m_uwnd_dvrtdx_pertremainder(LEVEL)
        #aa.f_m_vwnd_dvrtdy_pertremainder(LEVEL)
        pass

if PLOT:
    x1=aa.m_uwndbar_dvrtdxprime
    tcoord=x1.coord('time')
    time1=tcoord.units.num2date(tcoord.points[-1])
    time_constraint=iris.Constraint(time=time1)
    x2=x1.extract(time_constraint)

    qplt.contourf(x2)
    plt.gca().coastlines()
    
    plt.show()
