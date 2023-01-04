"""<Placeholder description. Replace after copying file.>

#######################################################################
This script was copied from the template plot_multipanel.py.

Always update the template plot_multipanel.py for greater
functionality, rather than creating multiple versions on different
branches.

If your resolution wavers, refer to xkcd: The General Problem.

"""

import copy
import datetime
import numpy as np
import os
import pickle
import string

import cartopy.crs as ccrs
import cftime
import iris
import iris.plot as iplt
import matplotlib as mpl
VIEW=False
if not VIEW:
    # Set mpl.use before importing matplotlib.pyplot
    # For production of final images set VIEW to False as interactive plotting
    # can change plot layout
    #mpl.use('Agg') # Non-interactive backends. 'Agg' for PNGs,'PDF','SVG','PS'
    pass
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
from panels import FigureSizeLocator
import pdb

import data_analysis as da
import info
import plotter

# Headers for use in calls to print
h1a='<<<=============================================================\n'
h1b='=============================================================>>>\n'
h2a='<<<---------------------------\n'
h2b='--------------------------->>>\n'

#==========================================================================

class Plot(object):

    """Plotting object.

    """

#==========================================================================

    def __init__(self,**descriptor):
        """Initialise.

        Set basic parameter values.
        """
        self.__dict__.update(descriptor)
        self.descriptor=descriptor

        self.ADHOC1=False # Should be False unless good reason. Search code for (ad hoc!) use.

        self.NROW=3; self.NCOL=1

        self.XVAR='longitude'; self.YVAR='latitude'
        self.XVAR1=0; self.XVAR2=360; self.YVAR1=-90; self.YVAR2=90

        #self.XVAR='longitude'; self.YVAR='Level'
        #self.XVAR1=0; self.XVAR2=360; self.YVAR1=50; self.YVAR2=1000

        #self.XVAR='time'; self.YVAR='Level'
        #self.XVAR1=cftime.DatetimeGregorian(2018,1,2,0,0);
        #self.XVAR2=cftime.DatetimeGregorian(2018,2,26,0,0);
        #self.YVAR1=10; self.YVAR2=1000
 
        #self.XVAR='time'; self.YVAR='latitude'
        #self.XVAR1=cftime.DatetimeGregorian(2014,1,15,0,0);
        #self.XVAR2=self.XVAR1+datetime.timedelta(days=15)
        #self.YVAR1=-11; self.YVAR2=-5

        #self.XVAR='longitude'; self.YVAR='time'
        #self.XVAR1=90-0.51; self.XVAR2=130+0.51
        #self.YVAR1=cftime._cftime.DatetimeGregorian(999,12,28,0,0);
        #self.YVAR2=cftime._cftime.DatetimeGregorian(1000,1,5,0,0);

        self.READ_DATA=True
        self.PLOT_COASTLINES='110m' # False, '110m', '50m', '10m'
        self.PROJTYPE='PlateCarree' # Only used for lat-lon. Options so far are 'PlateCarree', 'NorthPolarStereo'
        self.TIME_EXTRACT=True # If data for each panel is being extracted based on time (or lag) values
        self.INVERT_YAXIS=False
        self.YAXIS_LOG=False
        self.XLABELS_BOTTOM_ONLY=True
        self.YLABELS_LEFT_ONLY=True
        self.PLOT_LINES=False # False, 'all', or list of panel indices, e.g., [0,1,2,4,5]
        self.PLOT_MARKERS='all' # False, 'all', or list of panel indices, e.g., [0,1,2,4,5]
        self.PRINT_TITLE=False
        self.PANEL_LABELS=True
        self.PRINT_FILENAMES=True
        self.VIEW=VIEW
        self.SAVE_IMAGE=True
        self.DPI=False # Can be False or eg 50 (for quick look), 300 (manuscript)
        
        self.LABELS_FONTSIZE=11 # for tick and colorbar labels.  Try 10-28
        self.XTICKLABEL_ROTATION=0 # degrees

        # LOOPPANEL is an iterable that is looped over for each panel
        # Could be phases of MJO, level, time of day, etc.
        #self.LOOPPANEL=[-999] # dummy
        #self.LOOPPANEL=range(self.NROW*self.NCOL)
        #self.LOOPPANEL=['rmm006all'+str(xx) for xx in range(1,8+1)]
        self.LOOPPANEL=[0,10,20]
        #self.LOOPPANEL=['CCEK102E98-18-00UTC-and-M0002b','M0002b-and-not-CCEK102E98-18-00UTC']
        #self.LOOPPANEL=[500,]
        #self.LOOPPANEL=[ cftime.DatetimeGregorian(2018,9,13,19,30),
        #                 cftime.DatetimeGregorian(2018,9,14,1,30),
        #                 cftime.DatetimeGregorian(2018,9,14,7,30),
        #                 cftime.DatetimeGregorian(2018,9,14,13,30),
        #                 cftime.DatetimeGregorian(2018,9,14,19,30),
        #                 cftime.DatetimeGregorian(2018,9,15,1,30),
        #                 cftime.DatetimeGregorian(2018,9,15,7,30),
        #                 cftime.DatetimeGregorian(2018,9,15,13,30) ]
        #self.LOOPPANEL=list(range(0,21+1,3)) # hour of day
        #self.LOOPPANEL=['dvrtdt','source_dvrtdt','res_dvrtdt']
        #self.LOOPPANEL=list(range(1,8+1)) # MJO phases
        #t0=cftime.DatetimeGregorian(2011,11,21); tdel=datetime.timedelta(days=1)
        #self.LOOPPANEL=[t0+xx*tdel for xx in range(self.NROW*self.NCOL)]

        # Other variables for use in file names or data extraction
        #self.SEC_NO=1; self.LOC1=info.sections[self.SEC_NO]['location1']; self.LOC2=info.sections[self.SEC_NO]['location2']
        self.YEAR=1000; self.MONTH=1; self.DAY_START=1 # lagged composites day zero is 1 Jan 1000, by convention
        #self.YEAR=2018; self.MONTH=9#; self.DAY_START=11 # Typhoon Ompong
        #self.YEAR=2019; self.MONTH=3; self.DAY=1
        self.TDOMAINID='rmm001a-n2a5'

        self.IMAGEFILE=os.path.join(os.path.sep,'gpfs','home','e058','tmp','fig1.png')
        #self.IMAGEFILE=os.path.join(os.path.sep,'gpfs','home','e058','tmp',str(self.YEAR)+'-'+str(self.MONTH).zfill(2)+'-'+str(self.DAY).zfill(2)+'.png')
        #self.IMAGEFILE=os.path.join(os.path.sep,'gpfs','home','e058','tmp','fig_'+str(cftime.DatetimeGregorian(self.YEAR,self.MONTH,self.DAY_START))+'.png')
        #self.IMAGEFILE=os.path.join(os.path.sep,'gpfs','home','e058','tmp',self.LOC1+'-'+self.LOC2+'_'+str(self.XVAR1)+'.png')

        self.BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')
        #self.BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')
        #self.BASEDIR=os.path.join(os.path.sep,'gpfs','home','e058','home','data')
        
        #self.SUBDIR='std'
        self.SUBDIR='processed'

        # Use SINGLEFILE for one-off file names that do not fit data naming convention
        # Note, if file name relies on e.g., var_name or level then use
        # 'constant' attribute later, instead of SINGLEFILE
        self.SINGLEFILE1=False
        #self.SINGLEFILE1=os.path.join(self.BASEDIR,'gebco_sfc_ti','std','elevation_1.nc')
        self.SINGLEFILE2=False
        self.SINGLEFILE3=False
        self.SINGLEFILE4=False
        #self.SINGLEFILE4=os.path.join(self.BASEDIR,'ymc','bias_feb.nc')
        self.SINGLEFILE5=False
        self.SINGLEFILE6=False
        self.SINGLEFILE7=False

        # Set variables
        self.delta=1e-6

        # If this is run from __main__, descriptor should be an empty
        #   (dummy) dictionary
        # If this is called from another script, descriptor should contain
        #   information to set, e.g., the object(s) to be plotted. The code
        #   later should be modified to account for this on an individual 
        #   basis. 
        # The class attributes are updated here using descriptor, to allow
        #   for overwriting of attributes.
        # Example, descriptor could be used to overrid LOOPPANEL with cubes
        #   to be plotted (rather than read).
        self.__dict__.update(descriptor)
        self.descriptor=descriptor

        # Calculate derived variables
        self.npanel=self.NROW*self.NCOL
        #self.npanel=3
        #
        self.xvarrange=self.XVAR2-self.XVAR1
        self.yvarrange=self.YVAR2-self.YVAR1
        if self.XVAR=='longitude' and self.YVAR=='latitude':
            if self.PROJTYPE=='PlateCarree':
                self.panelratio=self.xvarrange/self.yvarrange
            elif self.PROJTYPE=='NorthPolarStereo':
                self.panelratio=1
            else:
                raise ValueError('Invalid PROJTYPE')
        else:
            #self.panelratio=1.5
            self.panelratio=8.4/10.8 # Baranowski et al. 2016a, Fig. 2
        #if 'time' in [self.XVAR,self.YVAR]:
        #    time1val=tcoord.points[0]
        #    time2val=tcoord.points[-1]
        #    print('time1val,time2val: {0!s}, {1!s}'.format(time1val,time2val))
        #    timerange=time2val-time1val
        #    if self.XVAR=='time':
        #        self.xvarrange=timerange
        #    if self.YVAR=='time':
        #        self.yvarrange=timerange

#==========================================================================

    def __repr__(self):
        return 'Plot'
    
#==========================================================================

    def __str__(self):
            ss=h1a+'Plot instance \n'+\
                'NROW, NCOL: {0.NROW!s}, {0.NCOL!s} \n'+\
                'XVAR, YVAR: {0.XVAR!s}, {0.YVAR!s} \n'+\
                'descriptor: {0.descriptor!s}, \n'
            ss+=h1b
            return ss.format(self)

#==========================================================================

    def set_variable_dictionaries(self):
        """For each variable, set list of parameters.

        factor can be False or a scaling number. NB consider convert_units instead of factor."""
        # Variable 1: to be contourf'd
        self.var1=dict([('plot',True),
                        ('source','ncepdoe_plev_d'),
                        ('var','zg'),
                        ('level',200),
                        ('colorbar',True),
                        ('filepre','_rac'),
                        ('factor',False),
                        ('constant',False),
                        ('basedir',self.BASEDIR),
                        ('subdir',self.SUBDIR),
                        ('singlefile',self.SINGLEFILE1)])
        # Variable 2: to be contour'd in thin contours
        self.var2=dict([('plot',False),
                        ('source','eqwaves'),
                        ('var','uwnd'),
                        ('level',1),
                        ('filepre',''),
                        ('factor',False),
                        ('constant',False),
                        ('basedir',self.BASEDIR),
                        ('subdir',self.SUBDIR),
                        ('singlefile',self.SINGLEFILE2)])
        # Variables 3a,3b: to be vector plotted
        self.var3a=dict([('plot',False),
                         ('source','era5mcw_sfc_h'),
                         ('var','uwnd'),
                         ('level',1),
                         ('filepre',''),
                         ('factor',False),
                         ('constant',False),
                         ('basedir',self.BASEDIR),
                         ('subdir',self.SUBDIR),
                         ('singlefile',self.SINGLEFILE3)])
        self.var3b=copy.copy(self.var3a)
        self.var3b['var']='vwnd'
        self.var3b['factor']=False # False if plotting u,v, e.g. -50 if v,omega
        # Variable 4: to be contour'd in thick contours
        self.var4=dict([('plot',False),
                        ('source','eqwaves'),
                        ('var','div'),
                        ('level',1),
                        ('filepre',''),
                        ('factor',1e6),
                        ('constant',False),
                        ('basedir',self.BASEDIR),
                        ('subdir',self.SUBDIR),
                        ('singlefile',self.SINGLEFILE4)])
        # Variables 5a,5b: to be vector plotted (second set of vectors)
        self.var5a=dict([('plot',False),
                         ('source','eqwaves'),
                         ('var','uwndchi'),
                         ('level',1),
                         ('filepre',''),
                         ('factor',False),
                         ('constant',True),
                         ('basedir',self.BASEDIR),
                         ('subdir',self.SUBDIR),
                         ('singlefile',self.SINGLEFILE5)])
        self.var5b=copy.copy(self.var5a)
        self.var5b['var']='vwndchi'
        self.var5b['factor']=False # False if plotting u,v, e.g. -50 if v,omega
        # Variables 6a,6b: to be vector plotted (third set of vectors)
        self.var6a=dict([('plot',False),
                         ('source','eqwaves'),
                         ('var','uwndpsi'),
                         ('level',1),
                         ('filepre',''),
                         ('factor',False),
                         ('constant',True),
                         ('basedir',self.BASEDIR),
                         ('subdir',self.SUBDIR),
                         ('singlefile',self.SINGLEFILE6)])
        self.var6b=copy.copy(self.var6a)
        self.var6b['var']='vwndpsi'
        self.var6b['factor']=False # False if plotting u,v, e.g. -50 if v,omega
        # Variable 7: to be contourf'd with transparency
        self.var7=dict([('plot',False),
                        ('source','imergmcw_sfc_30m'),
                        ('var','ppt'),
                        ('level',1),
                         ('colorbar',True),
                        ('filepre',''),
                        ('factor',False),
                        ('constant',False),
                        ('basedir',self.BASEDIR),
                        ('subdir',self.SUBDIR),
                        ('singlefile',self.SINGLEFILE7)])
        #
        self.dictlist=[self.var1,self.var2,self.var3a,self.var3b,self.var4,self.var5a,self.var5b,self.var6a,self.var6b,self.var7]

        # Set directories and initialise empty lists of cubes
        for dictc in self.dictlist:
            if dictc['plot']:
                dictc['dir']=os.path.join(dictc['basedir'],dictc['source'],dictc['subdir'])
                dictc['cubes']=[]
        print('# Completed set_variable_dictionaries.')

#==========================================================================

    def set_2d_arrays(self):
        """Set 2-d arrays."""

        self.panel_labels=[]
        for looppanelc in self.LOOPPANEL:
            for dictc in self.dictlist:
                if dictc['plot']:
                    if dictc['singlefile']:
                        file1=dictc['singlefile']
                    elif dictc['constant']:
                        #file1=os.path.join(dictc['dir'],dictc['var']+'_'+str(dictc['level'])+'_hovWKfiltEK_lat_-2.625_2.625_1998-01-01_2019-09-29.nc')
                        file1=os.path.join(dictc['dir'],dictc['var']+'_'+str(dictc['level'])+'_hov_lat_-2.625_2.625_1998-01-01_2019-09-29.nc')
                    else:
                        #file1=os.path.join(dictc['dir'],dictc['var']+'_'+str(dictc['level'])+dictc['filepre']+'_'+str(self.YEAR)+str(self.MONTH).zfill(2)+'.nc')
                        #file1=os.path.join(dictc['dir'],dictc['var']+'_'+str(dictc['level'])+dictc['filepre']+'_2018.nc') # from eg anncycle.py
                        file1=os.path.join(dictc['dir'],dictc['var']+'_'+str(dictc['level'])+dictc['filepre']+'_'+self.TDOMAINID+'_lag.nc') # from eg lagged_mean.py
                        #file1=os.path.join(dictc['dir'],dictc['var']+'_'+str(dictc['level'])+dictc['filepre']+'_'+str(looppanelc)+'_lag_ss_lat_-2.625_2.625.nc') # from eg lagged_mean.py
                        #file1=os.path.join(dictc['dir'],dictc['var']+'_'+str(dictc['level'])+'_CCEK102E01-20-00UTC-and-M0002a_lag.nc')
                        #file1=os.path.join(dictc['dir'],dictc['var']+'_'+str(dictc['level'])+'_rac_rmm001djf'+str(looppanelc)+'.nc')
                        #file1=os.path.join(dictc['dir'],dictc['var']+'_'+str(dictc['level'])+'_201401_section'+str(self.SEC_NO).zfill(3)+'.nc')
                        #file1=os.path.join(dictc['dir'],dictc['var']+'_'+str(dictc['level'])+'_'+looppanelc+'.nc')
                        #file1=os.path.join(dictc['dir'],dictc['var']+'_'+str(dictc['level'])+'_'+str(self.YEAR)+str(self.MONTH).zfill(2)+'.nc')
                        #file1=os.path.join(dictc['dir'],dictc['var']+'_'+str(dictc['level'])+'_'+str(self.YEAR)+str(self.MONTH).zfill(2)+'.nc')
                        #file1=os.path.join(dictc['dir'],dictc['var']+'_'+str(dictc['level'])+dictc['filepre']+'_'+str(looppanelc.year).zfill(4)+str(looppanelc.month).zfill(2)+'.nc')
                        #if dictc['var']=='ppt':
                        #    dummyc='mam98-16'
                        #else:
                        #    dummyc='mam79-17'
                        #file1=os.path.join(dictc['dir'],dictc['var']+'_'+str(dictc['level'])+'_'+dummyc+'.nc')
                    print('looppanelc,file1: {0!s}, {1!s}'.format(looppanelc,file1))
                    if self.READ_DATA:
                        print('Reading data')
                        x1=iris.load(file1)
                        d1=x1.concatenate_cube()
                    else:
                        # Put code here to set d1 to other data (not read from 
                        #   file), e.g., passed in descriptor.
                        #pass
                        print('Getting data from looppanelc')
                        d1=looppanelc
                    if dictc['factor']:
                        d1.data*=dictc['factor']
                    # Subset the input data for the range (XVAR1,XVAR2) and (YVAR1,YVAR2)
                    # Can use intersection and extract methods of cube to do this.
                    # intersection method is better (can't remember why)
                    # However, intersection method needs the relevant coordinate axis to have a modulus
                    # Longitude has a modulus (360.0). So does latitude for some data sets (360.0!!!)
                    # Time, Level do not, so use extract method with a constraint for these.
                    if (self.XVAR,self.YVAR)==('longitude','latitude'):
                        d2=d1.intersection(longitude=(self.XVAR1,self.XVAR2),latitude=(self.YVAR1,self.YVAR2))
                    elif (self.XVAR,self.YVAR)==('latitude','Level'):
                        d1a=d1.intersection(latitude=(self.XVAR1,self.XVAR2))
                        Level_constraint=iris.Constraint(Level=lambda cell: self.YVAR1 <=cell<= self.YVAR2)
                        d2=d1a.extract(Level_constraint)
                    elif (self.XVAR,self.YVAR)==('longitude','Level'):
                        d1a=d1.intersection(longitude=(self.XVAR1,self.XVAR2))
                        Level_constraint=iris.Constraint(Level=lambda cell: self.YVAR1 <=cell<= self.YVAR2)
                        d2=d1a.extract(Level_constraint)
                    elif (self.XVAR,self.YVAR)==('time','latitude'):
                        d1a=d1.intersection(latitude=(self.YVAR1,self.YVAR2))
                        time_constraint=iris.Constraint(time=lambda cell: self.XVAR1 <=cell<= self.XVAR2)
                        d2=d1a.extract(time_constraint)
                    elif (self.XVAR,self.YVAR)==('time','Level'):
                        Level_constraint=iris.Constraint(Level=lambda cell: self.YVAR1 <=cell<= self.YVAR2)
                        d1a=d1.extract(Level_constraint)
                        time_constraint=iris.Constraint(time=lambda cell: self.XVAR1 <=cell<= self.XVAR2)
                        d2=d1a.extract(time_constraint)
                    elif (self.XVAR,self.YVAR)==('longitude','time'):
                        d1a=d1.intersection(longitude=(self.XVAR1,self.XVAR2))
                        time_constraint=iris.Constraint(time=lambda cell: self.YVAR1 <=cell<= self.YVAR2)
                        pdb.set_trace()
                        d2=d1a.extract(time_constraint)
                    else:
                        raise ValueError('XVAR,YVAR combination not recognised')
                    #
                    if self.TIME_EXTRACT:
                        # Option to extract data based on times or lags
                        #timec=looppanelc
                        # timec=cftime.DatetimeGregorian(1,1,1,looppanelc)
                        #timec=cftime.DatetimeGregorian(self.YEAR,self.MONTH,self.DAY)
                        #timec=cftime.DatetimeGregorian(self.YEAR,self.MONTH,self.DAY_START,0)+looppanelc*datetime.timedelta(days=1)
                        timec=cftime.DatetimeGregorian(self.YEAR,self.MONTH,self.DAY_START,0)+looppanelc*datetime.timedelta(days=1)
                        #lagc=0
                        #timec=cftime.DatetimeGregorian(self.YEAR,self.MONTH,self.DAY_START,0)+lagc*datetime.timedelta(hours=1)
                        print('timec: {0!s}'.format(timec))
                        time_constraint=iris.Constraint(time=timec)
                        d2=d2.extract(time_constraint)
                    #
                    # Ad hoc further subsetting etc.
                    #d2=d2.extract(iris.Constraint(latitude=70))
                    #if d2.units=='degK': d2.convert_units('degC')
                    #if d2.units=='Pascal/s': d2.convert_units('hectoPascal day-1')
                    if d2.units=='mm hr-1': d2.convert_units('mm day-1')
                    #
                    dim_coord_names=[xx.name() for xx in d2.dim_coords]
                    if False and 'longitude' in [self.XVAR,self.YVAR]:
                        # Sections calculated using curved_section.py have time as one
                        # dim coord, and either latitude or longitude as the other dim
                        # coord and then the other one as an aux coord.  Both the lon
                        # and lat axes have the same length.  If necessary, swap them
                        # round here so the one that is to be used for plotting is the
                        # dim coord
                        if 'longitude' not in dim_coord_names:
                            print('# Promote longitude aux coord to dim coord at expense of latitude')
                            iris.util.promote_aux_coord_to_dim_coord(d2,'longitude')
                        elif 'latitude' not in dim_coord_names:
                            print('# Promote latitude aux coord to dim coord at expense of longitude')
                            iris.util.promote_aux_coord_to_dim_coord(d2,'latitude')
                    if 'time' in [self.XVAR,self.YVAR]:
                        print('Resetting time in time coordinate.')
                        plotter.cube_reset_time_coord(d2)
                    if d2.coords()[1].name()!=self.XVAR:
                        print('Transposing data coordinates.')
                        d2.transpose()
                    data_min=d2.data.min()
                    data_max=d2.data.max()
                    print('d2 min,max: {0!s}, {1!s}'.format(data_min,data_max))
                    unitsc=d2.units
                    print('units: {0!s}'.format(unitsc))
                    dictc['cubes'].append(d2)
                    #
            if self.TIME_EXTRACT:
                #self.panel_labels.append(str(timec)[10:16])
                #str1=str(timec)[8:10]+' ' +plotter.months_Jan[int(str(timec)[5:7])] +' ' +str(timec)[11:16]+' UTC'
                #str1=str(timec)[8:10]+' ' +plotter.months_Jan[int(str(timec)[5:7])]
                #self.panel_labels.append(str1)
                pass
        print('# Completed set_2d_arrays.')

#==========================================================================

    def set_ticks_and_labels(self):

        if 'longitude' in [self.XVAR,self.YVAR]:
            print('# Create longitude ticks and labels')
            label_int=30.
            self.lon_ticks=np.arange(-180,360+self.delta,label_int)
            self.lon_labels= [plotter.lonlat2string(xx,'lon',format=1) for xx in self.lon_ticks]
            print('lon_ticks: {0.lon_ticks!s}'.format(self))
            print('lon_labels: {0.lon_labels!s}'.format(self))
        
        if 'latitude' in [self.XVAR,self.YVAR]:
            print('# Create latitude ticks and labels')
            label_int=30.
            self.lat_ticks=np.arange(-90,90+self.delta,label_int)
            self.lat_labels= [plotter.lonlat2string(xx,'lat',format=1) for xx in self.lat_ticks]
            if self.ADHOC1:
                self.lat_ticks=[19.9 if val==20 else val for val in self.lat_ticks] # Overwrite as at edge of dataset
            print('lat_ticks: {0.lat_ticks!s}'.format(self))
            print('lat_labels: {0.lat_labels!s}'.format(self))

        if 'Level' in [self.XVAR,self.YVAR]:
            print('# Create Level ticks and labels')
            levels=[10,20,30,50,70,100,150,200,250,300,400,500,600,700,850,925,1000]
            self.level_ticks=np.array(levels)
            if self.YVAR=='Level' and self.YAXIS_LOG:
                self.level_ticks=np.log(self.level_ticks)
            self.level_labels=levels
            print('level_ticks: {0.level_ticks!s}'.format(self))
            print('level_labels: {0.level_labels!s}'.format(self))

        if 'time' in [self.XVAR,self.YVAR]:
            print('# Create time ticks and labels')
            tcoord=self.var1['cubes'][0].coord('time')
            tunits=tcoord.units
            dates=tunits.num2date(tcoord.points)
            self.time0=tunits.num2date(tcoord.points[0])
            #nn=5; self.time_ticks=[tunits.date2num(xx) for xx in dates if xx.hour==0 and divmod(xx.day-1,nn)[1]==0] # 1,n+1th,2n+1th etc of month (e.g., if n=5, 1st, 6th, 11th etc)
            nn=40; self.time_ticks=[tunits.date2num(xx) for xx in dates[::nn]] # every nn'th point
            #self.time_ticks=[tunits.date2num(xx) for xx in dates if xx.day==1 and xx.hour==0]
            self.time_labels=[str(tcoord.units.num2date(xx).day)+' '+plotter.months_Jan[tcoord.units.num2date(xx).month] for xx in self.time_ticks]
            #self.time_labels=[plotter.months_1_Jan[tcoord.units.num2date(xx).month] for xx in self.time_ticks]
            print('time_ticks: {0.time_ticks!s}'.format(self))
            print('time_labels: {0.time_labels!s}'.format(self))
        print('# Completed set_ticks_and_labels.')

#==========================================================================

    def create_overall_figure(self):
        """Create overall figure."""
        #loc=FigureSizeLocator(self.NROW,self.NCOL,figwidth=300,panelratio=self.panelratio,hsep=18,vsep=12,padleft=15,padright=10,padtop=10,padbottom=80,units='mm')
        self.loc=FigureSizeLocator(self.NROW,self.NCOL,figwidth=300,panelratio=self.panelratio,hsep=6,vsep=6,padleft=35,padright=10,padtop=10,padbottom=80,units='mm')
        self.fig=plt.figure(figsize=self.loc.figsize)
        print('# Completed create_overall_figure.')

#==========================================================================

    def set_projections(self):
        """Set up projections for latitude-longitude plots."""
        if self.PROJTYPE=='PlateCarree':
            self.proj0=ccrs.PlateCarree() # central_longitude is default of 0
            self.proj180=ccrs.PlateCarree(central_longitude=180)
            # Still not totally sure how the central_longitude argument is used.
            # Code below is from trial and error
            # Set projc as projection to be used for PlateCarree projections
            if self.XVAR=='longitude' and self.YVAR=='latitude':
                if 0<self.XVAR1<180 and 0<self.XVAR2<180: # Eastern hemisphere
                    self.projc=self.proj0
                elif 0<=self.XVAR1<180 and 180<self.XVAR2<=360: # Straddling dateline, with longitude in Western hemisphere going from 180 to 360, not -180 to 0
                    self.projc=self.proj180
                elif -180<=self.XVAR1<0 and 0<self.XVAR2<=180: # Straddling 0E, with longitude in Western hemisphere going from -180 to 0, not 180 to 360
                    self.projc=self.proj0
                else:
                    raise da.ToDoError('Find by trial and error how it works for other ranges!')
        elif self.PROJTYPE=='NorthPolarStereo':
            self.proj0=ccrs.NorthPolarStereo()
            self.projc=self.proj0
        else:
            raise ValueError('Invalid PROJTYPE')
        print('# Completed set_projections.')

#==========================================================================

    def set_plotting_properties(self):
        """Set contour levels etc."""
        if self.var1['plot']:
            print('# var1. Colour shading.')
            # Colour blind safe
            # Divergent colour: RdYlBu RdBu PuOr BrBG PiYG PRGn
            # Sequential single hue: Reds Blues Greens Greys Oranges Purples
            # Sequential multi hue: OrRd PuBu YlOrBr YlOrRd RdPu etc
            # Number of levels: 03-11 for divergent, 03-09 for sequential
            cmap1_name='brewer_'+'RdBu_11'
            self.cmap1=mpl_cm.get_cmap(cmap1_name)
            #self.cmap1=plotter.reverse_colormap(self.cmap1)
            #self.cmap1=mpl.colors.ListedColormap(['cyan','violet','blue','green','yellow','orange','red','pink'])
            self.constant_interval=True
            if self.constant_interval:
                cint1=10; ndigits1=1; clow1=-4.5*cint1; chigh1=-clow1
                #cint1=250; ndigits1=0; clow1=0; chigh1=clow1+9*cint1
                #cint1=2.5; ndigits1=1; clow1=2.5; chigh1=14*cint1+clow1
                self.levels1=plotter.levels_list(clow1,chigh1,cint1,zero=True,ndigits=ndigits1)
            else:
                # Non-constant contour interval requires a different colorbar set up later
                #self.levels1=[0.1,0.2,0.5,1,2,4,8,16,32,64]
                self.levels1=[1,2,4,8,16,32]
                self.cmap1.set_over('purple')
                self.cmap1.set_under('white')
                self.norm=mpl.colors.BoundaryNorm(self.levels1,self.cmap1.N)
            print('levels1: {0.levels1!s}'.format(self))
            self.var1['cint']=self.levels1[1]-self.levels1[0]

        if self.var2['plot']:
            print('# var2. Line contours (nominally thin).')
            cint2=0.25; ndigits2=2; chigh2=10*cint2; clow2=-chigh2
            self.levels2=plotter.levels_list(clow2,chigh2,cint2,zero=False,ndigits=ndigits2)
            #self.levels2=[0]
            #self.levels2=self.levels1
            print('levels2: {0.levels2!s}'.format(self))
            nlevels2=len(self.levels2)
            if divmod(nlevels2,2)[1]==0:
                # even
                nn=round(nlevels2/2)
                self.colors2=nn*('blue',)+nn*('red',)
                self.linestyles2=nn*('--',)+nn*('solid',)
                #self.linewidths2=nn*(1,)+nn*(1,)
            else:
                # odd
                nn=round((nlevels2-1)/2)
                self.colors2=nn*('blue',)+('black',)+nn*('red',)
                self.linestyles2=nn*('solid',)+('solid',)+nn*('solid',)
                self.linewidths2=nn*(1,)+(2,)+nn*(1,)
            self.colors2='black'
            #self.linestyles2='solid'
            self.linewidths2=1
            if len(self.levels2)>=2:
                self.var2['cint']=self.levels2[1]-self.levels2[0]

        if self.var3a['plot']:
            print('# var3. Vectors.')
            self.stridex=3
            self.stridey=3
            self.scale3=80 # scale3 is length of vector arrow in data units per scale_width (set later, usually mm). Bigger scale3 is, smaller vectors are.
            self.ref_arrow_size=5 # reference arrow for display (in units of input, e.g. m s-1). Independent of scale3
            self.color3='black' # 'black', 'pink', 'purple' are good
            self.var3a['cint']=self.var3b['cint']=self.ref_arrow_size

        if self.var4['plot']:
            print('# var4. Line contours (nominally thick).')
            cint4=0.5; ndigits4=1; chigh4=5*cint4; clow4=-chigh4
            self.levels4=plotter.levels_list(clow4,chigh4,cint4,zero=False,ndigits=ndigits4)
            #self.levels4=[-3,3]
            #self.levels4=self.levels2
            print('levels4: {0.levels4!s}'.format(self))
            nlevels4=len(self.levels4)
            if divmod(nlevels4,2)[1]==0:
                # even
                nn=round(nlevels4/2)
                self.colors4=nn*('blue',)+nn*('red',)
                self.linestyles4=nn*('--',)+nn*('solid',)
                #self.linewidths4=nn*(1,)+nn*(1,)
            else:
                # odd
                nn=round((nlevels4-1)/2)
                self.colors4=nn*('blue',)+('black',)+nn*('red',)
                self.linestyles4=nn*('solid',)+('solid',)+nn*('solid',)
                self.linewidths4=nn*(1,)+(2,)+nn*(1,)
            #self.colors4='blue'
            #self.linestyles4='solid'
            self.linewidths4=3
            self.var4['cint']=self.levels4[1]-self.levels4[0]

        if self.var5a['plot']:
            print('# var5. Vectors (second set).')
            # Assume plotting properties for second set of vectors are same as for first, except colour.
            self.color5='pink' # 'black', 'pink', 'purple' are good

        if self.var6a['plot']:
            print('# var6. Vectors (third set).')
            # Assume plotting properties for third set of vectors are same as for first, except colour.
            self.color6='purple' # 'black', 'pink', 'purple' are good

        if self.var7['plot']:
            print('# var7. Colour shading with transparency.')
            cmap7_name='brewer_'+'Blues_09'
            self.cmap7=mpl_cm.get_cmap(cmap7_name)
            self.constant_interval=True
            if self.constant_interval:
                cint7=2; ndigits7=0; clow7=0; chigh7=clow7+9*cint7
                self.levels7=plotter.levels_list(clow7,chigh7,cint7,zero=True,ndigits=ndigits7)
            else:
                # Non-constant contour interval requires a different colorbar set up later
                self.levels7=[1,2,4,8,16,32]
                self.cmap7.set_over('purple')
                self.cmap7.set_under('white')
                self.norm=mpl.colors.BoundaryNorm(self.levels7,self.cmap7.N)
            print('levels7: {0.levels7!s}'.format(self))
            self.var7['cint']=self.levels7[1]-self.levels7[0]
            self.alpha7=0.7 # Transparency for contourf

        print('# Completed set_plotting_properties.')

#==========================================================================

    def plot_figure(self):
        """Plot figure."""
        
        print('# Loop over panels')
        self.axlist=np.zeros((self.NROW,self.NCOL)).tolist()
        for ipanel in range(self.npanel):
            # subplot
            panel_order='down_first' # 'down_first' or 'right_first'
            if panel_order=='down_first':
                irow=divmod(ipanel,self.NROW)[1]
                icol=divmod(ipanel,self.NROW)[0]
            elif panel_order=='right_first':
                irow=divmod(ipanel,self.NCOL)[0]
                icol=divmod(ipanel,self.NCOL)[1]
            else:
                raise ValueError("panel_order must be either 'down_first' or 'right_first'")
            print('ipanel,irow,icol: {0!s}, {1!s} {2!s}'.format(ipanel,irow,icol))
            pos=self.loc.panel_position(irow,icol)
            if self.XVAR=='longitude' and self.YVAR=='latitude':
                axc=self.fig.add_axes(pos,projection=self.projc)
            else:
                axc=self.fig.add_axes(pos)
            self.axlist[irow][icol]=axc
            
            # plot axes ticks and determine tick labels
            axc.tick_params(direction='out')
            #axc.tick_params(length=10) # Override default tick length if labels overlap at corner
            axc.grid(True)
            if self.XVAR=='longitude' and self.YVAR=='latitude' and self.PLOT_COASTLINES:
                axc.coastlines(color='purple',linewidth=1,resolution=self.PLOT_COASTLINES)
            if self.XVAR=='longitude':
                if self.YVAR=='latitude':
                    axc.set_xticks(self.lon_ticks,crs=self.proj0)
                else:
                    axc.set_xticks(self.lon_ticks)
                xticklabels=self.lon_labels
            elif self.XVAR=='latitude':
                axc.set_xticks(self.lat_ticks)
                xticklabels=self.lat_labels
            elif self.XVAR=='time':
                axc.set_xticks(self.time_ticks)
                xticklabels=self.time_labels
            if self.YVAR=='latitude':
                if self.XVAR=='longitude':
                    axc.set_yticks(self.lat_ticks,crs=self.proj0)
                else:
                    axc.set_yticks(self.lat_ticks)
                yticklabels=self.lat_labels
            elif self.YVAR=='Level':
                axc.set_yticks(self.level_ticks)
                yticklabels=self.level_labels
            elif self.YVAR=='time':
                axc.set_yticks(self.time_ticks)
                yticklabels=self.time_labels
            # determine axes labels
            xlabel=ylabel=''
            #xlabel='Longitude'; ylabel='Pressure'
            #xlabel='Day ('+plotter.months_Jan[self.time0.month]+' '+str(self.time0.year)+'): tick mark at 0000 UTC'
            #xlabel='UTC '+str(icol*3)
            #ylabel=str(self.DAY_START+irow)+' '+plotter.months_Jan[self.MONTH]
            # plotting axes tick labels and labels
            ticklabels_fontsize=self.LABELS_FONTSIZE
            if not(self.XLABELS_BOTTOM_ONLY) or (self.XLABELS_BOTTOM_ONLY and irow==self.NROW-1) and not(self.PROJTYPE=='NorthPolarStereo'):
                axc.set_xticklabels(xticklabels,fontsize=ticklabels_fontsize,rotation=self.XTICKLABEL_ROTATION)
                axc.set_xlabel(xlabel,fontsize=ticklabels_fontsize)
            else:
                axc.set_xticklabels([])
                axc.tick_params(labelbottom='off')
            if not(self.YLABELS_LEFT_ONLY) or (self.YLABELS_LEFT_ONLY and icol==0) and not(self.PROJTYPE=='NorthPolarStereo'):
                axc.set_yticklabels(yticklabels,fontsize=ticklabels_fontsize)
                axc.set_ylabel(ylabel,fontsize=ticklabels_fontsize)
            else:
                axc.set_yticklabels([])
                axc.tick_params(labelleft='off')
            if self.INVERT_YAXIS:
                axc.invert_yaxis()

            # contourf variable 1
            if self.var1['plot']:
                cube1=self.var1['cubes'][ipanel]
                xmin=cube1.data.min(); xmax=cube1.data.max()
                print('cube1 min,max: {0!s}, {1!s}'.format(xmin,xmax))
                if self.XVAR=='longitude' and self.YVAR=='latitude':
                    cs1=iplt.contourf(cube1,levels=self.levels1,extend='both',cmap=self.cmap1)
                    # extend can be "neither", "both", "min", or "max"
                else:
                    # iplt.contourf does not produce time labels if time is a coord?
                    xdat1=cube1.coords()[1].points
                    ydat1=cube1.coords()[0].points
                    if self.YAXIS_LOG:
                        ydat1=np.log(ydat1)
                    cs1=axc.contourf(xdat1,ydat1,cube1.data,levels=self.levels1,extend='both',cmap=self.cmap1)

            # contourf variable 7
            if self.var7['plot']:
                cube7=self.var7['cubes'][ipanel]
                xmin=cube7.data.min(); xmax=cube7.data.max()
                print('cube7 min,max: {0!s}, {1!s}'.format(xmin,xmax))
                if self.XVAR=='longitude' and self.YVAR=='latitude':
                    cs7=iplt.contourf(cube7,levels=self.levels7,extend='max',cmap=self.cmap7,alpha=self.alpha7)
                    # extend can be "neither", "both", "min", or "max"
                else:
                    # iplt.contourf does not produce time labels if time is a coord?
                    xdat1=cube7.coords()[1].points
                    ydat1=cube7.coords()[0].points
                    if self.YAXIS_LOG:
                        ydat1=np.log(ydat1)
                    cs7=axc.contourf(xdat1,ydat1,cube7.data,levels=self.levels7,extend='both',cmap=self.cmap7,alpha=self.alpha7)

            # contour variable 4
            if self.var4['plot']:
                cube4=self.var4['cubes'][ipanel]
                xmin=cube4.data.min(); xmax=cube4.data.max()
                print('cube4 min,max: {0!s}, {1!s}'.format(xmin,xmax))
                levels4c,colors4c,linestyles4c,linewidths4c=plotter.contour_prune(cube4,self.levels4,self.colors4,self.linestyles4,self.linewidths4)
                if len(levels4c)>0:
                    if self.XVAR=='longitude' and self.YVAR=='latitude':
                        cs4=iplt.contour(cube4,levels=levels4c,colors=colors4c,linestyles=linestyles4c,linewidths=linewidths4c)
                    else:
                        xdat4=cube4.coords()[1].points
                        ydat4=cube4.coords()[0].points
                        if self.YAXIS_LOG:
                            ydat4=np.log(ydat4)
                        cs4=axc.contour(xdat4,ydat4,cube4.data,levels=levels4c,colors=colors4c,linestyles=linestyles4c,linewidths=linewidths4c)

            # contour variable 2
            if self.var2['plot']:
                cube2=self.var2['cubes'][ipanel]
                xmin=cube2.data.min(); xmax=cube2.data.max()
                print('cube2 min,max: {0!s}, {1!s}'.format(xmin,xmax))
                levels2c,colors2c,linestyles2c,linewidths2c=plotter.contour_prune(cube2,self.levels2,self.colors2,self.linestyles2,self.linewidths2)
                if len(levels2c)>0:
                    if self.XVAR=='longitude' and self.YVAR=='latitude':
                        cs2=iplt.contour(cube2,levels=levels2c,colors=colors2c,linestyles=linestyles2c,linewidths=linewidths2c)
                    else:
                        xdat2=cube2.coords()[1].points
                        ydat2=cube2.coords()[0].points
                        if self.YAXIS_LOG:
                            ydat2=np.log(ydat2)
                        cs2=axc.contour(xdat2,ydat2,cube2.data,levels=levels2c,colors=colors2c,linestyles=linestyles2c,linewidths=linewidths2c)

            # vector plot variable 3a,3b
            if self.var3a['plot']:
                var3_xdat=self.var3a['cubes'][ipanel].coord(self.XVAR).points
                var3_ydat=self.var3a['cubes'][ipanel].coord(self.YVAR).points
                if self.XVAR=='longitude' and self.YVAR=='latitude':
                    cs3=axc.quiver(var3_xdat[::self.stridex],var3_ydat[::self.stridey],self.var3a['cubes'][ipanel].data[::self.stridey,::self.stridex],self.var3b['cubes'][ipanel].data[::self.stridey,::self.stridex],color=self.color3,pivot='middle',transform=self.proj0,scale=self.scale3,scale_units='width')
                else:
                    cs3=axc.quiver(var3_xdat[::self.stridex],var3_ydat[::self.stridey],self.var3a['cubes'][ipanel].data[::self.stridey,::self.stridex],self.var3b['cubes'][ipanel].data[::self.stridey,::self.stridex],color=self.color3,pivot='middle',scale=self.scale3,scale_units='width')
                if ipanel==self.npanel-1:
                    # Plot reference arrow for last panel
                    left=axc.get_position().x0
                    right=axc.get_position().x1
                    ref_arrow_xpos=left+0.9*(right-left)
                    ref_arrow_yoffset=-0.65 # This value has to be surprisingly large (and negative)
                    ref_arrow_ypos=axc.get_position().y0+ref_arrow_yoffset
                    ref_arrow_text=str(self.ref_arrow_size)+' m s$^{-1}$'
                    #ref_arrow_text='[u]: '+str(self.ref_arrow_size)+' m s$^{-1}$. [$\omega$]: '+str(self.ref_arrow_size*86400/(self.var3b['factor']*-100))+' hPa day$^{-1}$'
                    ref_arrow_fontsize=self.LABELS_FONTSIZE
                    qk3=axc.quiverkey(cs3,ref_arrow_xpos,ref_arrow_ypos,self.ref_arrow_size,ref_arrow_text,transform=axc.transAxes,labelpos='W',fontproperties={'size':ref_arrow_fontsize})

            # vector plot variable 5a,5b
            if self.var5a['plot']:
                var5_xdat=self.var5a['cubes'][ipanel].coord(self.XVAR).points
                var5_ydat=self.var5a['cubes'][ipanel].coord(self.YVAR).points
                if self.XVAR=='longitude' and self.YVAR=='latitude':
                    cs5=axc.quiver(var5_xdat[::self.stridex],var5_ydat[::self.stridey],self.var5a['cubes'][ipanel].data[::self.stridey,::self.stridex],self.var5b['cubes'][ipanel].data[::self.stridey,::self.stridex],color=self.color5,pivot='middle',transform=self.proj0,scale=self.scale3,scale_units='width')
                else:
                    cs5=axc.quiver(var5_xdat[::self.stridex],var5_ydat[::self.stridey],self.var5a['cubes'][ipanel].data[::self.stridey,::self.stridex],self.var5b['cubes'][ipanel].data[::self.stridey,::self.stridex],color=self.color5,pivot='middle',scale=self.scale3,scale_units='width')

            # vector plot variable 6a,6b
            if self.var6a['plot']:
                var6_xdat=self.var6a['cubes'][ipanel].coord(self.XVAR).points
                var6_ydat=self.var6a['cubes'][ipanel].coord(self.YVAR).points
                if self.XVAR=='longitude' and self.YVAR=='latitude':
                    cs6=axc.quiver(var6_xdat[::self.stridex],var6_ydat[::self.stridey],self.var6a['cubes'][ipanel].data[::self.stridey,::self.stridex],self.var6b['cubes'][ipanel].data[::self.stridey,::self.stridex],color=self.color6,pivot='middle',transform=self.proj0,scale=self.scale3,scale_units='width')
                else:
                    cs6=axc.quiver(var6_xdat[::self.stridex],var6_ydat[::self.stridey],self.var6a['cubes'][ipanel].data[::self.stridey,::self.stridex],self.var6b['cubes'][ipanel].data[::self.stridey,::self.stridex],color=self.color6,pivot='middle',scale=self.scale3,scale_units='width')

            # Plot lines
            if self.PLOT_LINES and (self.PLOT_LINES=='all' or ipanel in self.PLOT_LINES):
                print('# Add lines')
                if True:
                    # Vertical line
                    x1=x2=110
                    tcoord=self.var1['cubes'][0].coord('time')
                    y1=tcoord.points[0]
                    y2=tcoord.points[-1]
                    if self.XVAR=='longitude' and self.YVAR=='latitude':
                        line0=mpl.lines.Line2D((x1,x2),(y1,y2),transform=self.proj0,linewidth=3,color='black')
                    else:
                        line0=mpl.lines.Line2D((x1,x2),(y1,y2),linewidth=3,color='black')
                    axc.add_line(line0)
                if False:
                    x1,y1=info.getlonlat('Christmas_Island')
                    x2,y2=info.getlonlat('Jakarta')
                    line1=mpl.lines.Line2D((x1,x2),(y1,y2),transform=self.proj0,linewidth=3,color='green')
                    axc.add_line(line1)
                if True:
                    # Read in CCKW database and extract trajectories that lie in 
                    # the current domain
                    self.file_traj_pickle=os.path.join(self.BASEDIR,self.var1['source'],'processed',self.var1['var']+'_'+str(self.var1['level'])+'_trajEK_lat_-2.625_2.625_1998-01-01_2019-09-29.pkl')
                    print('file_traj_pickle: {0.file_traj_pickle!s}'.format(self))
                    filec=open(self.file_traj_pickle,'rb')
                    self.trajectories=pickle.load(filec)
                    traj_in_domain=[]
                    ttunits=self.var1['cubes'][ipanel].coord('time').units
                    print('ttunits: {0!s}'.format(ttunits))
                    print('Trajectories that lie in domain:')
                    def testoverlap(x1,x2,y1,y2):
                        """Return true if interval [x1,x2] overlaps [y1,y2], otherwise return False."""
                        if x1<=y1 and x2>=y1 or x1>=y1 and x2<=y2 or x1<=y2 and x2>=y2:
                            return True
                        else:
                            return False
                    for keyc in self.trajectories.keys():
                        trajc=self.trajectories[keyc]
                        lon1=trajc['lons'][0]
                        lon2=trajc['lons'][-1]
                        time1=trajc['times'][0]
                        time2=trajc['times'][-1]
                        if testoverlap(lon1,lon2,self.XVAR1,self.XVAR2) and testoverlap(time1,time2,self.YVAR1,self.YVAR2):
                            print('keyc, lon1, lon2, time1, time2: {0!s}, {1!s}, {2!s}, {3!s}, {4!s}'.format(keyc,lon1,lon2,time1,time2))
                            traj_in_domain.append(trajc)
                            xx=[ trajc['lons'][ii] for ii in range(trajc['npts']) if self.XVAR1<=trajc['lons'][ii]<=self.XVAR2 and self.YVAR1<=trajc['times'][ii]<=self.YVAR2 ]
                            yy=[ ttunits.date2num(trajc['times'][ii]) for ii in range(trajc['npts']) if self.XVAR1<=trajc['lons'][ii]<=self.XVAR2 and self.YVAR1<=trajc['times'][ii]<=self.YVAR2 ]
                            linec=mpl.lines.Line2D(xx,yy,linewidth=3,color='red')
                            axc.add_line(linec)

            # Plot markers
            if self.PLOT_MARKERS and (self.PLOT_MARKERS=='all' or ipanel in self.PLOT_MARKERS):
                print('# Add markers')
                if False:
                    x1,y1=info.getlonlat('ELOX1')
                    axc.plot(x1,y1,color='purple',marker='o',markersize=10,transform=self.proj0)
                    x1,y1=info.getlonlat('ELOX2')
                    axc.plot(x1,y1,color='green',marker='d',markersize=10,transform=self.proj0)
                if True:
                    x1,y1=info.getlonlat('Padang')
                    axc.plot(x1,y1,color='green',marker='o',markersize=10,transform=self.proj0)

            # Set panel labels
            if self.PANEL_LABELS:
                self.panel_labels=['DJF','JJA'] # Can overwrite here
                #text=self.panel_labels[ipanel]
                #text=str(self.LOOPPANEL[ipanel])
                text='lag '+str(self.LOOPPANEL[ipanel])+' d'
                #text=str(self.LOOPPANEL[ipanel]); text=text.upper()[0]+text[1:]
                # Add (a), (b), etc to label;  chr(97) returns 'a', chr(98) returns 'b'
                text='('+chr(ipanel+97)+') '+text
                props=dict(boxstyle='round', facecolor='wheat', alpha=1.0)
                panel_label_position='tl' # 'bl', 'tl', 'br', 'tr'
                hoffset=0.05
                voffset=0.05
                if panel_label_position=='bl':
                    axc.text(hoffset,voffset,text,transform=axc.transAxes,fontsize=self.LABELS_FONTSIZE,horizontalalignment='left',verticalalignment='bottom',bbox=props)
                elif panel_label_position=='tl':
                    axc.text(hoffset,1-voffset,text,transform=axc.transAxes,fontsize=self.LABELS_FONTSIZE,horizontalalignment='left',verticalalignment='top',bbox=props)
                elif panel_label_position=='br':
                    axc.text(1-hoffset,voffset,text,transform=axc.transAxes,fontsize=self.LABELS_FONTSIZE,horizontalalignment='right',verticalalignment='bottom',bbox=props)
                elif panel_label_position=='tr':
                    axc.text(1-hoffset,1-voffset,text,transform=axc.transAxes,fontsize=self.LABELS_FONTSIZE,horizontalalignment='right',verticalalignment='top',bbox=props)

        # Set colorbar for var1
        if self.var1['plot'] and self.var1['colorbar']:
            print('# Colorbar for var1')
            colorbar_yoffset=0.05
            colorbar_position='horizontal_below_centred_figure'
            colorbar_fontsize=self.LABELS_FONTSIZE
            if colorbar_position=='horizontal_below_last_panel':
                axc=self.axlist[self.NROW-1][self.NCOL-1]
                bottom=axc.get_position().y0-colorbar_yoffset
                width=0.8*(axc.get_position().x1-axc.get_position().x0)
                middle=0.5*(axc.get_position().x1+axc.get_position().x0)
                left=middle-0.5*width
                height=0.015
            elif colorbar_position=='horizontal_below_centred_figure':
                axc=self.axlist[self.NROW-1][self.NCOL-1]
                bottom=axc.get_position().y0-colorbar_yoffset
                ax0=self.axlist[0][0]
                width=0.8*(axc.get_position().x1-ax0.get_position().x0)
                middle=0.5*(axc.get_position().x1+ax0.get_position().x0)
                left=middle-0.5*width
                height=0.015
            else:
                raise ValueError('Invalid colorbar_position.')
            cbaxes1=self.fig.add_axes([left,bottom,width,height])
            ndigits_cbar=1
            colorbar_labels='bottom' # 'bottom', 'top_and_bottom'
            if colorbar_labels=='bottom':
                # Lables on bottom of colorbar (default)
                xticks=[round(zz,ndigits_cbar) for zz in self.levels1] # all labels
            elif colorbar_labels=='top_and_bottom':
                # Alternate labels on top and bottom of colorbar
                ltop1=[round(zz,ndigits_cbar) for zz in self.levels1[::2]] # labels appear on top
                lbot1=[round(zz,ndigits_cbar) for zz in self.levels1[1:][::2]] # labels appear at bottom
                xticks=lbot1
                xmin=self.levels1[0]
                if cs1.extend in ['both','min']:
                    xmin-=self.var1['cint']
                xmax=self.levels1[-1]
                if cs1.extend in ['both','max']:
                    xmax+=self.var1['cint']
                for ii in ltop1:
                    strc=str(ii)
                    xfrac=(ii-xmin)/(xmax-xmin)
                    print(ii,xfrac,strc)
                    cbaxes1.text(xfrac, 1.3, strc, transform=cbaxes1.transAxes, va='bottom', ha='center',fontsize=colorbar_fontsize)
            if self.constant_interval:
                cbar1=plt.colorbar(cs1,cax=cbaxes1,orientation='horizontal',ticks=xticks,extendfrac='auto')
            else:
                sm=mpl.cm.ScalarMappable(cmap=self.cmap1,norm=self.norm)
                sm.set_array([]) # Not clear why this is needed
                cbar1=plt.colorbar(sm,cax=cbaxes1,orientation='horizontal',ticks=xticks,spacing='uniform',boundaries=[self.levels1[0]-self.var1['cint'],]+self.levels1+[self.levels1[-1]+self.var1['cint'],],extend='both',extendfrac='auto')
            cbar1.ax.tick_params(labelsize=colorbar_fontsize)
            #cbar1_label='Precipitation rate (mm day$^{-1}$)'
            #cbar1_label='Local Hadley circulation (kg m$^{-2}$ s$^{-1}$)'
            #cbar1_label='Vertical pressure velocity (Pa s$^{-1}$)'
            #cbar1_label='1000 hPa temperature anomaly ($^\circ$C)'
            #cbar1_label='Zonal wind (m s$^{-1}$)'
            #cbar1_label='Vorticity (x $10^{-6}$ s$^{-1}$)'
            #cbar1_label='Vorticity tendency (x $10^{-12}$ s$^{-2}$)'
            #cbar1_label='OLR (W m$^{-2}$)'
            cbar1_label='Geopotential height anomaly (m)'
            cbar1.set_label(cbar1_label,fontsize=colorbar_fontsize)

        # Set colorbar for var7
        if self.var7['plot'] and self.var7['colorbar']:
            print('# Colorbar for var7')
            colorbar_yoffset=0.07
            colorbar_position='horizontal_below_centred_figure'
            colorbar_fontsize=self.LABELS_FONTSIZE
            if colorbar_position=='horizontal_below_last_panel':
                axc=self.axlist[self.NROW-1][self.NCOL-1]
                bottom=axc.get_position().y0-colorbar_yoffset
                width=0.7*(axc.get_position().x1-axc.get_position().x0)
                middle=0.5*(axc.get_position().x1+axc.get_position().x0)
                left=middle-0.5*width
                height=0.015
            elif colorbar_position=='horizontal_below_centred_figure':
                axc=self.axlist[self.NROW-1][self.NCOL-1]
                bottom=axc.get_position().y0-colorbar_yoffset
                ax0=self.axlist[0][0]
                width=0.6*(axc.get_position().x1-ax0.get_position().x0)
                middle=0.5*(axc.get_position().x1+ax0.get_position().x0)
                left=middle-0.5*width
                height=0.015
            else:
                raise ValueError('Invalid colorbar_position.')
            cbaxes7=self.fig.add_axes([left,bottom,width,height])
            ndigits_cbar=1
            colorbar_labels='bottom' # 'bottom', 'top_and_bottom'
            if colorbar_labels=='bottom':
                # Lables on bottom of colorbar (default)
                xticks=[round(zz,ndigits_cbar) for zz in self.levels7] # all labels
            elif colorbar_labels=='top_and_bottom':
                # Alternate labels on top and bottom of colorbar
                ltop1=[round(zz,ndigits_cbar) for zz in self.levels1[::2]] # labels appear on top
                lbot1=[round(zz,ndigits_cbar) for zz in self.levels1[1:][::2]] # labels appear at bottom
                xticks=lbot1
                xmin=self.levels1[0]
                if cs7.extend in ['both','min']:
                    xmin-=self.var1['cint']
                xmax=self.levels1[-1]
                if cs7.extend in ['both','max']:
                    xmax+=self.var1['cint']
                for ii in ltop1:
                    strc=str(ii)
                    xfrac=(ii-xmin)/(xmax-xmin)
                    print(ii,xfrac,strc)
                    cbaxes7.text(xfrac, 1.3, strc, transform=cbaxes1.transAxes, va='bottom', ha='center',fontsize=colorbar_fontsize)
            if self.constant_interval:
                cbar7=plt.colorbar(cs7,cax=cbaxes7,orientation='horizontal',ticks=xticks,extendfrac='auto')
            else:
                sm=mpl.cm.ScalarMappable(cmap=self.cmap7,norm=self.norm)
                sm.set_array([]) # Not clear why this is needed
                cbar7=plt.colorbar(sm,cax=cbaxes7,orientation='horizontal',ticks=xticks,spacing='uniform',boundaries=[self.levels7[0]-self.var7['cint'],]+self.levels7+[self.levels7[-1]+self.var7['cint'],],extend='both',extendfrac='auto')
            cbar7.ax.tick_params(labelsize=colorbar_fontsize)
            cbar7_label='Precipitation rate (mm day$^{-1}$)'
            #cbar7_label='Local Hadley circulation (kg m$^{-2}$ s$^{-1}$)'
            #cbar7_label='Vertical pressure velocity (Pa s$^{-1}$)'
            #cbar7_label='1000 hPa temperature anomaly ($^\circ$C)'
            #cbar7_label='Zonal wind (m s$^{-1}$)'
            #cbar7_label='Vorticity (x $10^{-6}$ s$^{-1}$)'
            #cbar7_label='Vorticity tendency (x $10^{-12}$ s$^{-2}$)'
            #cbar7_label='OLR (W m$^{-2}$)'
            #cbar7_label='Elevation (m)'
            cbar7.set_label(cbar7_label,fontsize=colorbar_fontsize)

        # Set title
        if self.PRINT_TITLE:
            print('# Figure title')
            #fig_title='TRMM3B42 3-hourly precipitation: '+str(cftime.DatetimeGregorian(self.YEAR,self.MONTH,self.DAY_START))[:10]
            #fig_title='MAM mean: NCEP-DOE (u,$\omega$) averaged 10$^\circ$S-10$^\circ$N'
            #fig_title=plotter.months_January[MONTH]+' '+str(YEAR)
            #fig_title='GPM IMERG mean: '+self.TDOMAINID
            fig_title='Daily mean IMERG precipitation, ERA-5 surface wind: '+str(self.YEAR)+'-'+str(self.MONTH).zfill(2)+'-'+str(self.DAY).zfill(2)
            if False:
                fig_title='Anomalous: '# header text or ''
                for varc in self.dictlist:
                    if varc['plot']:
                        fig_title+='{0!s} {1!s}, {2!s} {3!s}; '.format(varc['var'],varc['level'],varc['cint'],varc['cubes'][0].units)
            print('fig_title: {0!s}'.format(fig_title))
            self.fig.text(0.5,0.99,fig_title,transform=self.fig.transFigure,fontsize=14,verticalalignment='top',horizontalalignment='center')

        if self.PRINT_FILENAMES:
            print('# Print name of calling script and image file')
            self.fig.text(0.99,0.998,os.path.basename(__file__)+' '+os.path.basename(self.IMAGEFILE),transform=self.fig.transFigure,fontsize=9,verticalalignment='top',horizontalalignment='right',color='gray')
        print('# Completed plot_figure.')

#==========================================================================

    def save_image(self):
        if self.SAVE_IMAGE:
            print('# Save image')
            print('IMAGEFILE {0!s}'.format(self.IMAGEFILE))
            if self.IMAGEFILE[-4:]=='.eps':
                print('Saving as eps with white space cropped.')
                self.fig.savefig(self.IMAGEFILE,bbox_inches='tight')
            elif self.DPI:
                print('Saving with dpi: {0!s}'.format(self.DPI))
                self.fig.savefig(self.IMAGEFILE,dpi=self.DPI)
            else:
                print('Saving with default dpi.')
                self.fig.savefig(self.IMAGEFILE)
        print('# Completed save_image.')

#==========================================================================

    def run(self):
        self.set_variable_dictionaries()
        self.set_2d_arrays()
        self.set_ticks_and_labels()
        self.create_overall_figure()
        self.set_projections()
        self.set_plotting_properties()
        self.plot_figure()
        if self.VIEW:
            print('# Plot')
            plt.show()
        self.save_image()

#==========================================================================

if __name__=='__main__':
    # Create instance of Plot object
    print('Running from __main__')
    descriptor={} # Empty dictionary for running from main
    aa=Plot(**descriptor)
    print(aa)
    
    # Run plotting methods
    aa.run()
