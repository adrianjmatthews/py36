"""Data analysis module using iris.

Hello world.

2018-04-25. New version for iris 2.

Classes that provide iris data i/o wrappers to data analysis methods,
often also using the iris module, for analysis of meteorological and
oceanographic gridded data.

Author: Adrian Matthews

Create documentation with pydoc -w data_analysis

Programming style:

Object oriented!

Class methods typically set an attribute(s) of the class instance.
They typically do not return an argument.

Printed output and information:

The __init__ method of each class has a 'verbose' keyword argument.

verbose=False (or 0) suppresses printed output

verbose=True (or 1) prints limited output (typically a statement that a
particular class attribute has been created)

verbose=2 prints extended output (typically the value of that class
attribute).

"""

from __future__ import division, print_function, with_statement # So can run this in python2

#if __name__ == '__main__':
# All import statements here, before class definitions
import copy
import datetime
import fnmatch
import os.path
import pdb
import pickle
import random
import shutil
import time

import cf_units
import cftime
import gsw # Gibbs TEOS-10 seawater routines
import iris
from iris.time import PartialDateTime
import math
import numpy as np
import scipy.signal # Note each module within the scipy package needs to be imported explicitly (import scipy# will not do it)
from windspharm.iris import VectorWind
from windspharm.tools import prep_data, recover_data
import spharm.spharm

import info
import mypaths

# Headers for use in calls to print
h1a='<<<=============================================================\n'
h1b='=============================================================>>>\n'
h2a='<<<---------------------------\n'
h2b='--------------------------->>>\n'

# var_name, standard_name/long_name pairs
# var_name is used for file names and as the variable name in netcdf files
# When iris extracts a cube from a netcdf file with a given name string, it
#   first looks for the standard_name, then long_name, then var_name.
# However, there is a prescribed list of allowed standard_name values
# 'tendency_of_atmosphere_relative_vorticity' (or any alternative) is not on
# the list, and an error will result if this is set as a standard_name.
# Use the iris.cube.Cube.rename(name) method to set the name of a cube
#   This sets the standard_name attribute if name is a valid standard_name
#   Otherwise it sets the long_name attribute to name.
# ie do not explicitly set the standard_name attribute as this will cause an
#   error if the desired name is not a valid standard_name
# rename method sets var_name to None if name is not a valid standard_name
# So var_name needs to be reset explicitly after rename has been called
# To set the name of cube1 to the name of cube 2, use
#   cube1.rename(cube2.name())
#   cube1.var_name=cube2.var_name
var_name2long_name={
    'bsiso1-1':'BSISO_1-1_index',
    'bsiso1-2':'BSISO_1-2_index',
    'bsiso1_amp':'BSISO_1_amplitude',
    'bsiso1_cat':'BSISO_1_category',
    'bsiso1_theta':'BSISO_1_phase_angle',
    'bsiso2-1':'BSISO_2-1_index',
    'bsiso2-2':'BSISO_2-2_index',
    'bsiso2_amp':'BSISO_2_amplitude',
    'bsiso2_cat':'BSISO_2_category',
    'bsiso2_theta':'BSISO_2_phase_angle',
    'chi':'atmosphere_horizontal_velocity_potential',
    'cphasex':'eastward_phase_speed',
    'cphasey':'northward_phase_speed',
    'cphasez':'upward_phase_speed',
    'div':'divergence_of_wind',
    'div_lambda':'divergence_of_wind_partitioned_eastward_component',
    'div_phi':'divergence_of_wind_partitioned_northward_component',
    'domegady_duwnddp':'meridional_derivative_of_lagrangian_tendency_of_air_pressure_times_pressure_derivative_of_zonal_wind',
    'dummy':'dummy_variable',
    'dvrtdt':'tendency_of_atmosphere_relative_vorticity',
    'duwnddx':'zonal_derivative_of_eastward_wind',
    'dvwnddx':'zonal_derivative_of_northward_wind',
    'dvwnddy':'meridional_derivative_of_northward_wind',
    'elevation':'height_above_reference_ellipsoid',
    'ew':'water_vapor_partial_pressure_in_air',
    'freq':'frequency',
    'ff':'coriolis_parameter',
    'harm':'integer_harmonic_in_time',
    'ke':'specific_kinetic_energy_of_air',
    'kk':'wavenumber_in_eastward_direction',
    'lat':'latitude',
    'lhfd':'surface_downward_latent_heat_flux',
    'ltt':'linear_temporal_trend',
    'll':'wavenumber_in_northward_direction',
    'lon':'longitude',
    'lwrd':'surface_net_downward_longwave_flux',
    'mf_lambda':'upward_mass_flux_of_air_partitioned_eastward_component',
    'mf_phi':'upward_mass_flux_of_air_partitioned_northward_component',
    'mm':'wavenumber_in_upward_direction',
    'mslp':'air_pressure_at_sea_level',
    'mltt':'ocean_mixed_layer_thickness_defined_by_temperature',
    'mu':'atmosphere_overturning_potential',
    'm_beta_vwnd':'minus_meridional_derivative_of_coriolis_parameter_times_northward_wind',
    'm_domegadx_dvwnddp':'minus_zonal_derivative_of_lagrangian_tendency_of_air_pressure_times_pressure_derivative_of_northward_wind',
    'm_duwnddy':'minus_meridional_derivative_of_eastward_wind',
    'm_ff_div':'minus_coriolis_parameter_times_atmosphere_relative_vorticity',
    'm_omega_dvrtdp':'minus_lagrangian_tendency_of_air_pressure_times_pressure_derivative_of_atmosphere_relative_vorticity',
    'm_uwnd_dvrtdx':'minus_eastward_wind_times_zonal_derivative_of_atmosphere_relative_vorticity',
    'm_uwndbar_dvrtdxbar':'minus_time_mean_eastward_wind_times_time_mean_zonal_derivative_of_atmosphere_relative_vorticity',
    'm_uwndbar_dvrtdxprime':'minus_time_mean_eastward_wind_times_time_perturbation_zonal_derivative_of_atmosphere_relative_vorticity',
    'm_uwndprime_dvrtdxbar':'minus_time_perturbation_eastward_wind_times_time_mean_zonal_derivative_of_atmosphere_relative_vorticity',
    'm_uwndprime_dvrtdxprime':'minus_time_perturbation_eastward_wind_times_time_perturbation_zonal_derivative_of_atmosphere_relative_vorticity',
    'm_vwnd_dvrtdy':'minus_northward_wind_times_meridional_derivative_of_atmosphere_relative_vorticity',
    'm_vwndbar_dvrtdybar':'minus_time_mean_northward_wind_times_time_mean_meridional_derivative_of_atmosphere_relative_vorticity',
    'm_vwndbar_dvrtdyprime':'minus_time_mean_northward_wind_times_time_perturbation_meridional_derivative_of_atmosphere_relative_vorticity',
    'm_vwndprime_dvrtdybar':'minus_time_perturbation_northward_wind_times_time_mean_meridional_derivative_of_atmosphere_relative_vorticity',
    'm_vwndprime_dvrtdyprime':'minus_time_perturbation_northward_wind_times_time_perturbation_meridional_derivative_of_atmosphere_relative_vorticity',
    'm_vrt_div':'minus_divergence_of_wind_times_atmosphere_relative_vorticity',
    'm_vrtbar_divbar':'minus_time_mean_divergence_of_wind_times_time_mean_atmosphere_relative_vorticity',
    'm_vrtbar_divprime':'minus_time_perturbation_divergence_of_wind_times_time_mean_atmosphere_relative_vorticity',
    'm_vrtprime_divbar':'minus_time_mean_divergence_of_wind_times_time_perturbation_atmosphere_relative_vorticity',
    'm_vrtprime_divprime':'minus_time_perturbation_divergence_of_wind_times_time_perturbation_atmosphere_relative_vorticity',
    'nhfd':'surface_downward_heat_flux_in_sea_water',
    'olr':'toa_outgoing_longwave_flux',
    'omega':'lagrangian_tendency_of_air_pressure',
    'omega_lambda':'lagrangian_tendency_of_air_pressure_partitioned_eastward_component',
    'omega_phi':'lagrangian_tendency_of_air_pressure_partitioned_northward_component',
    'omegaf':'angular_frequency',
    'ppt':'lwe_precipitation_rate',
    'pa':'air_pressure',
    'phi':'geopotential',
    'psfc':'surface_air_pressure',
    'psi':'atmosphere_horizontal_streamfunction',
    'psioc_lambda':'eastward_component_of_vector_atmosphere_overturning_streamfunction',
    'psioc_phi':'northward_component_of_vector_atmosphere_overturning_streamfunction',
    'pv':'ertel_potential_vorticity',
    'res_dvrtdt':'residual_tendency_of_atmosphere_relative_vorticity',
    'rhum':'relative_humidity',
    'rmm1':'RMM_1_index',
    'rmm2':'RMM_2_index',
    'rmm_amp':'RMM_amplitude',
    'rmm_cat':'RMM_category',
    'sa':'sea_water_absolute_salinity',
    'shfd':'surface_downward_sensible_heat_flux',
    'shum':'specific_humidity',
    'source_dvrtdt':'total_source_of_tendency_of_atmosphere_relative_vorticity',
    'ss':'integer_zonal_wavenumber',
    'ssft':'sea_surface_foundation_temperature',
    'sst':'sea_surface_temperature',
    'swp':'sea_water_pressure',
    'swpd':'sea_water_potential_density',
    'swrd':'surface_net_downward_shortwave_flux',
    'ta':'air_temperature',
    'taux':'surface_downward_eastward_stress',
    'tauy':'surface_downward_northward_stress',
    'theta':'air_potential_temperature',
    'tsc':'sea_water_conservative_temperature',
    'uwnd':'eastward_wind',
    'uwndchi':'irrotational_component_of_eastward_wind',
    'uwndpsi':'nondivergent_component_of_eastward_wind',
    'vrt':'atmosphere_relative_vorticity',
    'vrt_horiz_adv':'minus_eastward_wind_times_zonal_derivative_of_atmosphere_relative_vorticity_minus_northward_wind_times_meridional_derivative_of_atmosphere_relative_vorticity',
    'vrt_stretch':'minus_divergence_of_wind_times_atmospheric_total_vorticity',
    'vrt_tilt':'minus_zonal_derivative_of_lagrangian_tendency_of_air_pressure_times_pressure_derivative_of_northward_wind_plus_meridional_derivative_of_lagrangian_tendency_of_air_pressure_times_pressure_derivative_of_zonal_wind',
    'vwnd':'northward_wind',
    'vwndchi':'irrotational_component_of_northward_wind',
    'vwndpsi':'nondivergent_component_of_northward_wind',
    'vwndptap':'product_of_northward_wind_anomaly_and_air_temperature_anomaly',
    'wndspd':'wind_speed',
    'wnddir':'wind_from_direction',
    'wwnd':'upward_air_velocity',
    'zg':'geopotential_height',
    }

# Create inverse dictionary
def invert_dict(d):
    return dict([(v,k) for k,v in d.items()])
long_name2var_name=invert_dict(var_name2long_name)

#==========================================================================

def source_info(aa):
    """Create attributes of an object based on the source attribute.

    aa is an object.

    Attributes created are:

    data_source:  e.g., 'ncepdoe', 'olrinterp'
    level_type:  e.g., 'plev' for pressure level, 'toa' for top of atmosphere
    frequency: e.g., 'd' for daily, '3' for 3-hourly

    outfile_frequency: e.g., 'year' or 'month'
    wildcard: e.g., '????' or '??????'
    timedelta: datetime.timedelta object corresponding to the frequency attribute
    calendar: 'gregorian' or '360_day'
    
    This function is typically called in the initialisation of class instances in this module.
    """
    # Split source attribute string using underscores as separators
    xx=aa.source.split('_')
    if len(xx)!=3:
        raise RuntimeError("source attribute '{0.source!s}' must have three parts separated by underscores".format(aa))
    aa.data_source=xx[0]
    aa.level_type=xx[1]
    aa.frequency=xx[2]
    # Check data_source attribute is valid
    valid_data_sources=['era5trp','era5plp','era5bar','era5mcw','erainterim','erainterimEK1','erainterimNEK1','erainterimNEK1T42','erainterimEK2','erainterimEK3','imergplp','imergmcw','imergmts','imergmt2','imergnpl','imergnp2','imergtrm','imergtrmp1','ncepdoe','ncepdoegg','ncepncar','olrcdr','olrinterp','ostial4nrttrp','ostial4reptrp','sg579m031oi01','sg534m031oi01','sg532m031oi01','sg620m031oi01','sg613m031oi01','sgallm031oi01','sstrey','trmm3b42v7','trmm3b42v7p1','trmm3b42v7p2','trmm3b42v7p3','trmm3b42v7p4','tropflux','hadgem2esajhog']
    if aa.data_source not in valid_data_sources:
        raise UserWarning('data_source {0.data_source!s} not valid'.format(aa))
    # Set outfile_frequency attribute depending on source information
    if aa.source in ['erainterim_sfc_d','erainterim_sfc_6h','erainterim_plev_6h','erainterimEK1_plev_6h','erainterimNEK1_plev_6h','erainterimNEK1T42_plev_6h','erainterimEK2_plev_6h','erainterimEK3_plev_6h','erainterim_plev_d','ncepdoe_plev_6h','ncepdoe_plev_d','ncepdoe_sfc_d','ncepdoegg_zlev_d','ncepdoe_zlev_d','ncepncar_plev_d','ncepncar_sfc_d','olrcdr_toa_d','olrinterp_toa_d','sstrey_sfc_7d','sg579m031oi01_zlev_h','sg534m031oi01_zlev_h','sg532m031oi01_zlev_h','sg620m031oi01_zlev_h','sg613m031oi01_zlev_h','sgallm031oi01_zlev_h','sstrey_sfc_d','tropflux_sfc_d','hadgem2esajhog_plev_d']:
        aa.outfile_frequency='year'
        aa.wildcard='????'
    elif aa.source in ['imergplp_sfc_30m','imergmcw_sfc_30m','imergmcw_sfc_dt','imergmts_sfc_30m','imergmt2_sfc_30m','imergnpl_sfc_30m','imergnp2_sfc_30m','imergtrm_sfc_30m','imergtrm_sfc_3h','imergtrmp1_sfc_3h','trmm3b42v7_sfc_3h','trmm3b42v7p1_sfc_3h','trmm3b42v7p2_sfc_3h','trmm3b42v7_sfc_d','trmm3b42v7p1_sfc_d','trmm3b42v7p3_sfc_d','trmm3b42v7p4_sfc_d','era5trp_plev_h','era5plp_plev_h','era5plp_sfc_h','era5bar_sfc_h','era5mcw_sfc_h','era5mcw_sfc_d','ostial4nrttrp_sfc_d','ostial4reptrp_sfc_d']:
        aa.outfile_frequency='month'
        aa.wildcard='??????'
    else:
        raise UserWarning('Need to specify outfile_frequency for this data source.')
    # timedelta attribute
    if aa.frequency[-1]=='m':
        if aa.frequency=='m':
            aa.timedelta=datetime.timedelta(minutes=1)
        else:
            aa.timedelta=datetime.timedelta(minutes=int(aa.frequency[:-1]))
    elif aa.frequency[-1]=='h':
        if aa.frequency=='h':
            aa.timedelta=datetime.timedelta(hours=1)
        else:
            aa.timedelta=datetime.timedelta(hours=int(aa.frequency[:-1]))
    elif aa.frequency[-1]=='d':
        if aa.frequency=='d':
            aa.timedelta=datetime.timedelta(days=1)
        else:
            aa.timedelta=datetime.timedelta(days=int(aa.frequency[:-1]))
    else:
        raise ToDoError('Need to code up for different frequency attribute.')

    # Set calendar attribute depending on data_source
    if aa.data_source in ['hadgem2esajhog']:
        aa.calendar='360_day'
    else:
        aa.calendar='gregorian'
            
    # Printed output
    if aa.verbose:
        ss=h2a+'source_info.  Created attributes: \n'+\
            'data source: {0.data_source!s} \n'+\
            'level_type: {0.level_type!s} \n'+\
            'frequency: {0.frequency!s} \n'+\
            'outfile_frequency: {0.outfile_frequency!s} \n'+\
            'calendar: {0.calendar!s} \n'+\
            'timedelta: {0.timedelta!s} \n'+\
            'wildcard: {0.wildcard!s} \n'+h2b
        print(ss.format(aa))
                
#==========================================================================

def clean_callback(cube,field,filename):
    """Deletes many attributes on iris load.
    
    Problem.  iris concatenate and merge (to create a single cube from
    a cube list) is very picky and will fail if there are any
    mismatching metadata between the cubes.  This function removes
    attributes from the time coordinate and basic metadata that
    typically fall foul of this.  These attributes are not useful
    anyway.

    Data providers (such as NOAA) sometimes change the attributes
    (trivially) in the middle of their data sets, so these all need to
    be stripped out.

    Ideally, this clean callback will only be used in preprocess.py,
    and all subsequent analysis will not need it.  Hence attributes
    that are built up during analysis will be self-consistent and will
    not need removing.
    
    Usage: as an argument in iris load  (...,callback=clean_callback).
    
    """
    # Delete the problem attributes from the time coordinate:
    # 'actual_range' expected to be different between cubes, so delete
    # 'coordinate_defines' usually set to 'start' or 'point'
    #     'start' indicates that e.g., the 0000 UTC time stamp refers to the
    #         start of the day over which the daily average (eg NCEP data)
    #         has been carried out over.  Unfortunately, NCEP erroneously
    #         change this to 'point' in 2015.  Need to delete it.
    for coordc in cube.coords():
        if coordc.name() in ['time','t']:
            time_coord=coordc
    att_list=['actual_range','coordinate_defines']
    for attribute in att_list:
        if attribute in time_coord.attributes:
            del cube.coord('time').attributes[attribute]
    # Or set the attributes dictionary of the time coordinate to empty:
    #cube.coord('time').attributes = {}
    # 
    # Similarly for latitude and longitude coordinaes
    att_list2=['_ChunkSizes','_CoordinateAxisType','comment','valid_max','valid_min']
    for coordc in cube.coords():
        if coordc.name() in ['latitude','longitude']:
            for attribute in att_list2:
                if attribute in coordc.attributes:
                    del coordc.attributes[attribute]
    
    # Similarly delete  some of the main attributes
    att_list=['actual_range','history','unpacked_valid_range','references',
              'References','dataset_title','title','long_name',
              'Conventions','GRIB_id','GRIB_name','comments','institution',
              'least_significant_digit','level_desc','parent_stat',
              'platform','precision','source','statistic','title','var_desc',
              'NCO','creation_date','invalid_units','metodology',
              'producer_agency','time_range','website','short_name',
              'description','dataset','_NCProperties','valid_max','valid_min',
              'tracking_id','table_id','date','time','name',
              'FROM_ORIGINAL_FILE__Metadata_Conventions',
              'FROM_ORIGINAL_FILE__geospatial_lat_resolution',
              'FROM_ORIGINAL_FILE__geospatial_lat_units',
              'FROM_ORIGINAL_FILE__geospatial_lon_resolution',
              'FROM_ORIGINAL_FILE__geospatial_lon_units',
              'FROM_ORIGINAL_FILE__netcdf_version_id',
              'FROM_ORIGINAL_FILE__platform',
              'FROM_ORIGINAL_FILE__product_version',
              'FROM_ORIGINAL_FILE__northernmost_latitude',
              'FROM_ORIGINAL_FILE__southernmost_latitude',
              'FROM_ORIGINAL_FILE__westernmost_longitude',
              'FROM_ORIGINAL_FILE__easternmost_longitude',
              'Metadata_Conventions',
              'geospatial_lat_resolution',
              'geospatial_lat_units',
              'geospatial_lon_resolution',
              'geospatial_lon_units',
              'geospatial_lat_max',
              'geospatial_lat_min',
              'geospatial_lon_max',
              'geospatial_lon_min',
              'netcdf_version_id',
              'platform',
              'product_version',
              'northernmost_latitude',
              'southernmost_latitude',
              'westernmost_longitude',
              'easternmost_longitude',
              'History',
              '_ChunkSizes',
              '_CoordSysBuilder',
              'acknowledgment',
              'cdm_data_type',
              'comment',
              'creator_email',
              'creator_name',
              'creator_url',
              'date_created',
              'file_quality_level',
              'gds_version_id',
              'id',
              'keywords',
              'keywords_vocabulary',
              'license',
              'metadata_link',
              'naming_authority',
              'processing_level',
              'project',
              'publisher_email',
              'publisher_name',
              'publisher_url',
              'reference',
              'sensor',
              'spatial_resolution',
              'standard_name_vocabulary',
              'start_time',
              'stop_time',
              'summary',
              'time_coverage_end',
              'time_coverage_start',
              'uuid',]
    for attribute in att_list:
        if attribute in cube.attributes:
            del cube.attributes[attribute]

    # A cube also has some "attributes" that are not in the attributes dictionary
    # Set these to empty strings
    cube.long_name=''

    # Set cell methods to empty tuple
    cube.cell_methods=()

#==========================================================================

def compare_cube_attributes(cube1,cube2):
    """Compare attributes of two cubes and find any non matches. """

    att1=cube1.attributes
    att2=cube2.attributes
    def compare_attributes(att1,att2):
        flag=True
        for keyc in att1:
            if keyc not in att2:
                flag=False
                raise UserWarning(str(keyc)+'is not an attribute of cube2')
            elif att2[keyc]!=att1[keyc]:
                flag=False
                raise UserWarning('Attributes not identical: '+str(keyc)+': '+str(att2[keyc])+' != '+str(att1[keyc]))
    print('Checking attributes in cube1')
    flag1=compare_attributes(att1,att2)
    print('Checking attributes in cube2')
    flag2=compare_attributes(att2,att1)
    print(flag1,flag2)
    return flag1 and flag2

#==========================================================================

def create_cube(array,oldcube,new_axis=False,new_var_name=False):
    """Create an iris cube from a numpy array and attributes from an old cube.

    Inputs:

    <array> must be a numpy array.

    <oldcube>  must be an iris cube.

    <new_axis> is either False, or an iris cube axis (e.g., a time axis).

    <new_var_name> is either False, or a string variable name, e.g., 'psi'

    Output:

    <newcube> is an iris cube of the same shape as <array>

    Usage:

    If <new_axis> is False then array and oldcube must have the same
    shape and dimensions.  <newcube> is then simply an iris cube
    created from <array> and the coordinate axes and attributes of
    <oldcube>.

    If <new_axis> is an iris cube axis (e.g., a time axis), then
    <array> and <oldcube> must have the same number of dimensions, but
    one of the axes (the one that will correspond to <new_axis>) can
    be a different length.  For example, <array> could have shape
    (30,73,144), <oldcube> could have shape (365,73,144) with
    dimensions of time,latitude,longitude, and <new_axis> would be a
    time axis of length 30.  <newcube> would then be a cube with data
    from <array>, the time axis of <new_axis> and latitude and
    longitude axes from <oldcube>, plus all the attributes of
    <oldcube>.

    If <new_var_name> is False, the new cube has the same name as the
    old cube, otherwise it is renamed.
    
    """
    if new_axis:
        # new_axis is an iris coordinate axis
        # Find the index of the corresponding axis in oldcube
        # (assumes the new axis is a dim_coord not a aux_coord)
        coord_name=new_axis.standard_name
        print('coord_name: {0!s}'.format(coord_name))
        coord_names=[dimc.standard_name for dimc in oldcube.dim_coords]
        print('coord_names: {0!s}'.format(coord_names))
        if coord_name not in coord_names:
            raise UserWarning('coord_name is not in coord_names')
        coord_index=coord_names.index(coord_name)
        print('coord_index: {0!s}'.format(coord_index))


        # Create a list of two-lists, each of form [coord,index] for
        # dim_coord's and aux_coord's
        # Note that aux_coord's do not have to be assigned to a dimension
        #   in which case they have an empty tuple for index
        dim_coord_names=[xx.name() for xx in oldcube.dim_coords]
        aux_coord_names=[xx.name() for xx in oldcube.aux_coords]
        print('dim_coord_names: {0!s}'.format(dim_coord_names))
        print('aux_coord_names: {0!s}'.format(aux_coord_names))
        dim_coords_and_dims=[ [oldcube.coord(xx),oldcube.coord_dims(xx)] for xx in dim_coord_names ]
        aux_coords_and_dims=[ [oldcube.coord(xx),oldcube.coord_dims(xx)] for xx in aux_coord_names ]
        # Overwrite the coord_index'th axis with new_axis two-list
        dim_coords_and_dims[coord_index]=[new_axis,coord_index]
    else:
        # new_axis is False.
        # Check that array and oldcube have same dimensions
        if array.shape!=oldcube.shape:
            raise UserWarning('Shape of array and oldcube must match.')
        #
        # Create a list of two-lists, each of form [coord,index] for
        # dim_coord's and aux_coord's
        # Note that aux_coord's do not have to be assigned to a dimension
        #   in which case they have an empty tuple for index
        dim_coord_names=[xx.name() for xx in oldcube.dim_coords]
        aux_coord_names=[xx.name() for xx in oldcube.aux_coords]
        print('dim_coord_names: {0!s}'.format(dim_coord_names))
        print('aux_coord_names: {0!s}'.format(aux_coord_names))
        dim_coords_and_dims=[ [oldcube.coord(xx),oldcube.coord_dims(xx)] for xx in dim_coord_names ]
        aux_coords_and_dims=[ [oldcube.coord(xx),oldcube.coord_dims(xx)] for xx in aux_coord_names ]
    # Create cube
    newcube=iris.cube.Cube(array,units=oldcube.units,attributes=oldcube.attributes,cell_methods=oldcube.cell_methods,dim_coords_and_dims=dim_coords_and_dims,aux_coords_and_dims=aux_coords_and_dims)
    # Name cube
    if new_var_name:
        new_long_name=var_name2long_name[new_var_name]
        newcube.rename(new_long_name)
        newcube.var_name=new_var_name
    else:
        newcube.rename(oldcube.name())
        newcube.var_name=oldcube.var_name

    return newcube

#==========================================================================

def conv_float32(old_array):
    """Convert a numpy array type to float32 ('single precision').

    Default numpy array type is float64 (64 bit, or 'double
    precision').  This leads to a doubling of disk space needed to
    store data, compared to float32 (32 bit, 'single precision'), with
    only spurious increase in precision.

    It is not a problem having float64 type numbers in memory, only
    when they are written to disk.  Hence, just call this function as
    a final wrapper for a numpy array that is about to be converted to
    an iris cube.  Do not bother with all the previous intermediate
    steps, as float64 is insidious and can easily creep back into to
    your arrays somewhere.

    Input: old_array should be a numpy array.

    Returns new_array, a numpy array of type float32.

    """

    new_array=old_array.astype(type('float32',(np.float32,),{}))
    return new_array

#==========================================================================

def block_times(aa,verbose=False):
    """Set start and end times for current block.

    Input <aa> is an object with the following attributes:

    outfile_frequency:  e.g., 'year' or 'month'
    year: integer for current year, e.g., 2016
    month: integer for current month, in range 1 to 12.
    calendar

    Returns time,time2 which are datetime.datetime objects (if
    aa.calendar is 'gregorian') or cftime.Datetime360Day objects (if
    aa.calendar is '360_day') for the start and end times for the
    current block.  These are then typically used to create an iris
    constraint on time for extracting data.

    If outfile_frequency is 'year', start and end times are 00 UTC on
    1 Jan of current year, and 1 second before 00 UTC on 1 Jan the
    following year.

    If outfile_frequency is 'month', start and end times are 00 UTC on
    1st of current month and year, and 1 second before 00 UTC on 1st
    of the following month.

    """
    timedelta_second=datetime.timedelta(seconds=1)
    #timedelta_day=datetime.timedelta(days=1)
    if aa.calendar not in ['gregorian','360_day']:
        raise UserWarning('calendar not valid.')
    if aa.outfile_frequency=='year':
        if aa.calendar=='gregorian':
            time1=datetime.datetime(aa.year,1,1) # 00:00:00 on 1 Jan
            #time2=datetime.datetime(aa.year+1,1,1)-timedelta_second # 23:59:59 on 31 Dec
            time2=datetime.datetime(aa.year,12,31,23,59,59) # 23:59:59 on 31 Dec
        elif aa.calendar=='360_day':
            time1=cftime.Datetime360Day(aa.year,1,1) # 00:00:00 on 1 Jan
            time2=cftime.Datetime360Day(aa.year,12,30,23,59,59) # 23:59:59 on 30 Dec
    elif aa.outfile_frequency=='month':
        if aa.calendar=='gregorian':
            time1=datetime.datetime(aa.year,aa.month,1) # 00:00:00 on 1st of month
            if aa.month!=12:
                time2=datetime.datetime(aa.year,aa.month+1,1)-timedelta_second # 23:59:59 on last of month
            else:
                time2=datetime.datetime(aa.year+1,1,1)-timedelta_second  # 23:59:59 on last of month
        elif aa.calendar=='360_day':
            time1=cftime.Datetime360Day(aa.year,aa.month,1) # 00:00:00 on 1st of month
            time2=cftime.Datetime360Day(aa.year,aa.month,30,23,59,59) # 23:59:59 on 30th of month
    else:
        raise ToDoError("Need to write code for other outfile_frequency.")
    if verbose:
        ss=h2a+'block_times\n'+\
            ' calendar: {0!s}\n time1: {1!s}\n time2: {2!s}\n'+h2b
        print(ss.format(aa.calendar,time1,time2))
    return time1,time2

#==========================================================================

def replace_wildcard_with_time(aa,string1):
    """Replace wild card question marks with current time.

    Input <aa> is an object with the following attributes:

    outfile_frequency:  e.g., 'year' or 'month'
    wildcard: e.g., '????' or '??????'
    year: integer for current year, e.g., 2016
    month: integer for current month, in range 1 to 12.

    Input <string1> is a string, usually a filename, which contains
    wild card question marks, e.g., ???? or ??????.

    The wild card question marks are replaced with the current date,
    e.g., YYYY or YYYYMM, depending on aa.outfile_frequency.

    Output is <string2>, typically a filename referring to a specific
    date or range of dates (e.g., year and month).

    """
    # Set wild card string
    if aa.outfile_frequency=='year':
        x2=str(aa.year)
    elif aa.outfile_frequency=='month':
        x2=str(aa.year)+str(aa.month).zfill(2)
    else:
        raise ToDoError('Need to code up for different outfile_frequency')
    # Replace wild card characters
    string2=string1.replace(aa.wildcard,x2)

    return string2

#==========================================================================

def iter_generator(input):
    """Convert <input> to an iterator and return it.

    If <input> is an integer, convert it to a list so it can be
    iterated over.  Otherwise, return unchanged.

    Typically, YEAR or MONTH is passed to this function.  They may
    already be in the form of a range-type iterator or a list, or just a single
    integer value.  This function converts the single integer value to
    an iterator (a list).  """
    if type(input)==type(1979):
        # input is an integer.  Convert to a list
        output=[input]
    else:
        # Return input unchanged
        output=input
    return output

#==========================================================================

def lat_direction(cube_in,direction,verbose=True):
    """Make cube latitude run north to south or south to north.

    Inputs:

    <cube_in> is an iris cube

    <direction> is a string, either 's2n' or 'n2s'

    Use iris.util.reverse if needed.

    Return <cube_out>, whose latitude axis runs according to
    <direction>.
    """
    # Check direction is valid
    if direction not in ['s2n','n2s']:
        raise UserWarning('Direction is invalid.')
    # Find direction of latitude coordinate of input cube
    lat_coord=cube_in.coord('latitude')
    lat_beg=lat_coord.points[0]
    lat_end=lat_coord.points[-1]
    if lat_beg<lat_end:
        input_direction='s2n'
    else:
        input_direction='n2s'
    # If input cube does not have desired latitude direction, reverse it
    if direction!=input_direction:
        print('Changing latitude direction.')
        # Find index of latitude coordinate
        dim_coord_names=[xx.var_name for xx in cube_in.dim_coords]
        lat_index=dim_coord_names.index('latitude')
        # Reverse latitude
        cube_out=iris.util.reverse(cube_in,lat_index)
        # Check direction of output cube
        lat_coord=cube_out.coord('latitude')
        lat_beg=lat_coord.points[0]
        lat_end=lat_coord.points[-1]
        if lat_beg<lat_end:
            output_direction='s2n'
        else:
            output_direction='n2s'
        if direction!=output_direction:
            raise UserWarning('Ouput direction does not match direction.')
        # Printed output
        if verbose:
            ss=h2a+'name: {4!s} \n'+\
                'lat_direction. \n'+\
                'direction, input_direction: {0!s}, {1!s} \n'+\
                'lat_index: {2!s} \n'+\
                'direction, output_direction: {0!s}, {3!s} \n'+h2b
            print(ss.format(direction,input_direction,lat_index,output_direction,cube_in.name()))
    else:
        print('Not changing latitude direction.')
        cube_out=cube_in

    return cube_out
    
#==========================================================================

def truncate(cube_in,truncation):
    """Spectrally truncate (triangular truncation) iris cube.

    Inputs:

    <cube_in> is an iris cube.

    <truncation> is an integer, to triangular truncate to in spectral
    space.

    Outputs:

    Return <cube_out>, the spectrally truncated cube.
    """

    # Find value of south2north
    south2north=f_south2north(cube_in,verbose=True)
    # Create spoof uwnd and vwnd cubes to initiate VectorWind instance.
    uwnd=create_cube(cube_in.data,cube_in,new_var_name='uwnd')
    vwnd=create_cube(cube_in.data,cube_in,new_var_name='vwnd')
    # Create VectorWind instance
    ww=VectorWind(uwnd,vwnd)
    # Truncate input cube
    cube_out=ww.truncate(cube_in,truncation=truncation)
    # Add cell method to describe truncation
    cm=iris.coords.CellMethod('point',coords=['latitude','longitude'],comments='truncated to T'+str(truncation))
    cube_out.add_cell_method(cm)
    if south2north:
        # Output of windspharm method.  Latitude runs north to south
        # Reverse latitude direction (to be consistent with input cube)
        cube_out=lat_direction(cube_out,'s2n')
    return cube_out

#==========================================================================

def inverse_laplacian(cube_in):
    """Calculate inverse Laplacian of an iris cube.

    Inputs:

    <cube_in> is an iris cube.

    Outputs:

    Return <cube_out>, the inverse Laplacian.
    """

    # Find value of south2north
    south2north=f_south2north(cube_in,verbose=True)
    # Create spoof uwnd and vwnd cubes to initiate VectorWind instance.
    uwnd=create_cube(cube_in.data,cube_in,new_var_name='uwnd')
    vwnd=create_cube(cube_in.data,cube_in,new_var_name='vwnd')
    # Create VectorWind instance
    ww=VectorWind(uwnd,vwnd)
    # Get the Spharmt instance
    S = ww._api.s  # (just vw.s if you used windspharm.standard.VectorWind)
    # Transform your scalar field to spectral domain, you can't use metadata here
    scalar_field, shape_info = prep_data(cube_in.data, "tyx")
    scalar_field_spec = S.grdtospec(scalar_field)
    # Compute the inverse Laplacian
    ilap_spec = spharm.spharm._spherepack.invlap(scalar_field_spec, S.rsphere)
    ilap = recover_data(S.spectogrd(ilap_spec), shape_info)
    ilap_cube = uwnd.copy()
    ilap_cube.data = ilap
    ilap_cube.rename("inverse_laplacian")
    print(ilap_cube)
    cube_out=ilap_cube
    if south2north:
        # Output of windspharm method.  Latitude runs north to south
        # Reverse latitude direction (to be consistent with input cube)
        cube_out=lat_direction(cube_out,'s2n')
    return cube_out

#==========================================================================

def gradient(cube_in):
    """Calculate horizontal gradient of an iris cube.

    Inputs:

    <cube_in> is an iris cube.

    Outputs:

    Return <cube_gradx>,<cube_grady>, the eastward and northward gradient 
    fields.
    """
    # Find value of south2north
    south2north=f_south2north(cube_in,verbose=True)
    # Create spoof uwnd and vwnd cubes to initiate VectorWind instance.
    uwnd=create_cube(cube_in.data,cube_in,new_var_name='uwnd')
    vwnd=create_cube(cube_in.data,cube_in,new_var_name='vwnd')
    # Create VectorWind instance
    ww=VectorWind(uwnd,vwnd)
    # Calculate eastward and northward gradients of input cube
    cube_gradx,cube_grady=ww.gradient(cube_in)
    if south2north:
        # Output of windspharm method.  Latitude runs north to south
        # Reverse latitude direction (to be consistent with input cube)
        cube_gradx=lat_direction(cube_gradx,'s2n')
        cube_grady=lat_direction(cube_grady,'s2n')
    return cube_gradx,cube_grady

#==========================================================================

def concatenate_cube(cubelist):
    """concatenate_cube cube list when one cube has singleton time dimension.

    Input is <cubelist>, an iris CubeList

    Output is <cube>, a concatenate_cube'd version of <cubelist>

    The problem: If a cube list contains a cube that has a singleton
    time coordinate, then iris will have relegated that time
    coordinate to a scalar coordinate when it extracted the cube.  The
    other cubes will have time coordinates as dimension coordinates.
    If the cubelist concatenate_cube method is used, it will fail because
    of this mismatch of dimensions.

    Example:
    x1=iris.load('uwnd_200_????.nc') # Load data stored in single file per year
    time1=datetime.datetime(1994,12,31)
    time2=datetime.datetime(1995,1,6)
    time_constraint=iris.Constraint(time=lambda cell: time1<=cell<=time2)
    x3=x2.concatenate_cube()

    x2 is a cube list containing two cubes.  The first cube only has
    one time (31 Dec 1994) and the time coordinate will have been
    relegated to a scalar coordinate by the extract method.  The
    second cube has 6 time values (1-6 Jan 1995) and its time
    coordinate will be a dimension coordinate.  The
    x2.concatenate_cube() will then fail.

    Solution: Wherever this is likely to happen, replace the standard
    cube list concatenate_cube method, i.e., in the example above,
    replace

    x3=x2.concatenate_cube() 

    with a call to this function:

    x3=concatenate_cube(x2)

    This function uses iris.util.new_axis() to promote the singleton
    scalar time coordinate to a dimension coordinate.  The cubes can
    then be concatenate_cube'd using the standard cube list
    concatenate_cube method.

    """
    # If there is only one cube in the cube list, concatenate_cube as usual
    if len(cubelist)==1:
        cube=cubelist.concatenate_cube()
        return cube
    # If there is more than one cube in the cube list, check for singleton
    #   time dimensions
    for index in range(len(cubelist)):
        cubec=cubelist[index]
        time_coord=cubec.coord('time')
        if time_coord.shape==(1,):
            print('concatenate_cube.  Promoting singleton time dimension.')
            cubec=iris.util.new_axis(cubec,'time')
            cubelist[index]=cubec
    # Concatenate into single cube
    cube=cubelist.concatenate_cube()

    return cube

#==========================================================================

def cube_section(cube_in,xdat,ydat,verbose=False):
    """Calculate an irregular section through a cube.

    Called from e.g., curved_section.py

    Inputs:

    cube_in: input iris cube

    xdat,ydat : 2-tuples.  First element is string of coordinate name.
    Second element is list of coordinate values.  Both coordinate
    names must be coordinates of cube_in.

    e.g., xdat is ('longitude',[4.3, 5.7, 6.1])
          ydat is ('latitude', [13.1, 12.6, 5.6])

    Outputs:

    cube_out: iris cube, with data from cube_in interpolated onto the
    xdat,ydat values.

    """
    xdat_coord=xdat[0]
    xdat_vals=xdat[1]
    ydat_coord=ydat[0]
    ydat_vals=ydat[1]
    if verbose:
        ss=h1a+'section \n'+\
            'xdat_coord: {0!s}\n'+\
            'xdat_vals: {1!s}\n'+\
            'ydat_coord: {2!s}\n'+\
            'ydat_vals: {3!s}\n'+h1b
        print(ss.format(xdat_coord,xdat_vals,ydat_coord,ydat_vals))
    if len(xdat_vals)!=len(ydat_vals):
        raise UserWarning('xdat_vals and ydat_vals must have same number of values.')
    x1=iris.cube.CubeList([])
    for ival in range(len(xdat_vals)):
        xdat_valc=xdat_vals[ival]
        ydat_valc=ydat_vals[ival]
        samples=[(xdat_coord,xdat_valc),(ydat_coord,ydat_valc)]
        x2=cube_in.interpolate(samples,iris.analysis.Linear())
        x1.append(x2)
    cube_out=x1.merge_cube()

    return cube_out

#==========================================================================

def create_counter_from_mask(cube_in,verbose=False):
    """Create counter array based on data mask.

    Input:

    cube_in: input iris cube.

    Outputs:

    counter: numpy array with same shape as cube_in.  Filled with 1s
    and 0s: 1 if the relevant value of cube_in exists, 0 if it does
    not (i.e., is masked).

    cube_out: output iris cube with same shape as cube_in. Filled with
    mixture of data values and 0s: data value if the relevant value of
    cube_in exists, 0 if it does not (i.e., is masked).

    """
    try:
        mask=cube_in.data.mask
        if mask.any():
            message="Mask exists"
            counter=np.where(mask,0,1)
            x1=np.where(mask,0,cube_in.data)
            cube_out=create_cube(conv_float32(x1),cube_in)
        else:
            message="Mask exists, but is False! Counter is all 1's"
            counter=np.ones(cube_in.data.shape)
            cube_out=cube_in
    except AttributeError:
        message="No mask. Counter is all 1's"
        counter=np.ones(cube_in.data.shape)
        cube_out=cube_in
    print('create_counter_from_mask:',message)
    print('counter.shape, sum: {0!s}, {1!s}'.format(counter.shape,np.sum(counter)))

    return counter,cube_out

#==========================================================================

def find_npd(source):
    """Find npd (number per day) from source string.

    Output : npd
    """
    xx=source.split('_')
    source_frequency=xx[2]
    if source_frequency[-1]=='d':
        if source_frequency=='d':
            npd=1
        else:
            raise UserWarning('Cannot have frequency less than one per day.')
    elif source_frequency[-1]=='h':
        if source_frequency=='h':
            npd=24
        else:
            npd=int(24/int(source_frequency[:-1]))
    elif source_frequency[-1]=='m':
        if source_frequency=='m':
            npd=1440
        else:
            npd=int(1440/int(source_frequency[:-1]))
    else:
        raise ToDoError('Code up for source frequencies not in minutes or hours.')
    print('source_frequency,npd: {0!s}, {1!s}'.format(source_frequency,npd))
    return npd

#==========================================================================

def set_time_constraint(time1,time2,calendar='gregorian',verbose=False):
    """Set iris time constraint from datetime.datetime like object(s) according to calendar. 
    
    Inputs:
    <time1> a datetime.datetime or cftime.Datetime360Day object
    <time2> a datetime.datetime or cftime.Datetime360Day object or False
    <calendar> string to specify calendar

    Always call this function to set an iris time constraint, rather
    than explicitly setting iris.Constraint in the code.

    Return an iris.Constraint object: time constraint between time1
    and time2 inclusive, or just equal to time1 if time2 is False.

    This function is used as a single call to allow an iris time
    constraint to be set for data on either a gregorian or 360_day
    calendar. Data with a gregorian calendar needs datetime.datetime
    objects to set the time constraint. Date with a 360_day calendar
    needs cftime.Datetime360Day objects to set the time
    constraint. This function can take datetime.datetime objects and
    then internally converts them to cftime.Datetime360day objects if
    the calendar is 360_day before creating the iris time
    constraint. Alternatively, if the calendar is 360_day then can just pass cftime.Datetime360day objects to this function and it will create a
    """
    if calendar=='gregorian':
        # Check inputs are datetime.datetime objects
        # Or cftime._cftime.real_datetime object, which seem to have crept in somewhere
        # and appear to be equivalent to datetime.datetime objects?
        dummydatetime=datetime.datetime(1900,1,1)
        dummyrealdatetime=cftime._cftime.real_datetime(1900,1,1)
        if type(time1) not in [type(dummydatetime),type(dummyrealdatetime)] or (time2 and type(time2)!=type(dummydatetime)):
            raise UserWarning('calendar is gregorian. Inputs must be datetime.datetime or cftime._cftime.real_datetime objects.')
    elif calendar=='360_day':
        # Create cftime.Datetime360Day objects to set time constraint
        # NB inputs can be datetime.datetime or already cftime.Datetime360Day objects
        time1a=create_Datetime360Day(time1)
        time2a=create_Datetime360Day(time2)
    else:
        raise UserWarning('calendar not recognised.')
    if time2:
        # Set time constraint to be between time1 and time2, inclusive
        if calendar=='gregorian':
            time_constraint=iris.Constraint(time=lambda cell: time1<=cell<=time2)
        elif calendar=='360_day':
            time_constraint=iris.Constraint(time=lambda cell: time1a<=cell.point<=time2a)
    else:
        # Set time constraint to be equal to time1 +/- small increment to
        # allow for rounding errors in time coordinates
        timedelta1=datetime.timedelta(seconds=1)
        if calendar=='gregorian':
            time_constraint=iris.Constraint(time=lambda cell: time1-timedelta1<=cell<=time1+timedelta1)
        if calendar=='360_day':
            time_constraint=iris.Constraint(time=lambda cell: time1a-timedelta1<=cell<=time1a+timedelta1)
    if verbose:
        ss=h2a+'set_time_constraint \n'+\
            ' calendar: {0!s}\n time1: {1!s}\n time2: {2!s}\n'+h2b
        print(ss.format(calendar,time1,time2))
    return time_constraint

#==========================================================================

def create_Datetime360Day(time1):
    """Return a cftime.Datetime360Day object from a datetime like object.
    
    Input: <time1> a datetime like object, e.g., a datetime.datetime
    or cftime.Datetime360Day object, or False.
    """
    if time1:
        return cftime.Datetime360Day(time1.year,time1.month,time1.day,time1.hour,time1.minute,time1.second)
    else:
        return time1

#==========================================================================

def create_datetime(time1):
    """Return a datetime.datetime object from a datetime like object.
    
    Input: <time1> a datetime like object, e.g., a datetime.datetime
    or cftime.Datetime360Day object, or False.
    """
    if time1:
        return datetime.datetime(time1.year,time1.month,time1.day,time1.hour,time1.minute,time1.second)
    else:
        return time1

#==========================================================================

def change_time_stamp_from_12_to_00(self,verbose=False):
    """Time stamp is at 12 UTC. Change to 00 UTC by subtracting 0.5 (days)."""
    tcoord=self.cube.coord('time')
    time_units=tcoord.units
    if 'day' not in time_units.name:
        raise UserWarning('Expecting time units in days.')
    if divmod(tcoord.points[0],1)[1]!=0.5:
        raise UserWarning('Expecting time stamp to be at 12 UTC.')
    # Cannot simply subtract 0.5 because of round off, leading to later errors
    # Use divmod to achieve same end
    #time_val=tcoord.points-0.5
    time_val=divmod(tcoord.points,1)[0]
    tcoord2=iris.coords.DimCoord(time_val,standard_name='time',units=time_units)
    if self.verbose==2:
        print('tcoord: {0!s}'.format(tcoord))
        print('tcoord2: {0!s}'.format(tcoord2))
    self.cube.remove_coord('time')
    self.cube.add_dim_coord(tcoord2,0)
    print('Subtracted 0.5 days from time coord so time is now at 00 UTC not 12 UTC.')

#==========================================================================

def f_longitude_average(cube_in,nave):
    """Return longitude average ("regrid") of iris cube.
    
    Used to change (reduce) longitude resolution of data, for
    subsequent data analysis.
    
    Input:
    
    <cube_in> : iris cube that has longitude as one dimension.
    
    <nave> : integer number of longitudes to average over. E.g.,
    if nave is 4 and the input cube has a longitude resolution of
    1/4 degree, the output will have a longitude resolution of 1
    degree, each new longitude data point being the average of
    nave source data points.
    
    Output:
    
    <cube_out> : iris cube with data averaged (regridded) over
    longitude.
    
    Method: Based on method in f_time_average(method=2). Slice the
    input cube numpy array over longitude with steps of nave, to
    get nave separated cubes. Add them together and divide by
    nave. This is fast.  Depends on longitude values being equally
    spaced.
    """
    # Check length of longitude coord is a multiple of nave
    lon1_coord=cube_in.coord('longitude')
    nlon1=len(lon1_coord.points)
    if divmod(nlon1,nave)[1]!=0:
        raise UserWarning('nlon1 must be divisible by nave.')
    nlon2=divmod(nlon1,nave)[0]
    print('nlon1,nlon2,nave: {0!s},{1!s},{2!s}'.format(nlon1,nlon2,nave))
    # Check input data is well formed
    dim_coord_names=[xx.var_name for xx in cube_in.dim_coords]
    lon1_index=dim_coord_names.index('longitude')
    print('lon1_index: {0!s}'.format(lon1_index))
    # Slice numpy array into nave slices by longitude.
    # Then calculate longitude means
    x1=cube_in.data
    print('x1.shape: {0!s}'.format(x1.shape))
    if lon1_index+1==len(x1.data.shape):
        # This code requires longitude to be the last index
        x2=[x1[...,ii:nlon1:nave] for ii in range(nave)] 
    else:
        raise ToDoError('Use a modification of line above if longitude is not last dimension.')
    x3=x2[0]
    print('x3.shape: {0!s}'.format(x3.shape))
    for xx in x2[1:]:
        x3=x3+xx
        print('xx.shape: {0!s}'.format(xx.shape))
    x4=x3/nave
    print('x4.shape: {0!s}'.format(x4.shape))
    # Create new longitude axis for longitude averaged data
    # Calculate new longitudes by same method as calcuated longitude
    # averaged data
    lon1_vals=lon1_coord.points
    if lon1_index+1==len(x1.data.shape):
        # This code requires longitude to be the last index
        y2=[lon1_vals[...,ii:nlon1:nave] for ii in range(nave)] 
    else:
        raise ToDoError('Use a modification of line above if longitude is not last dimension.')
    y3=y2[0]
    for yy in y2[1:]:
        y3=y3+yy
    lon2_vals=y3/nave
    print('lon1_vals: {0!s}'.format(lon1_vals))
    print('lon2_vals: {0!s}'.format(lon2_vals))
    lon2_coord=iris.coords.DimCoord(lon2_vals,standard_name='longitude',units=lon1_coord.units)
    # Create new cube of longitude averaged data
    x11=create_cube(x4,cube_in,new_axis=lon2_coord)
    cm=iris.coords.CellMethod('point','longitude',comments='regridded in longitude by averaging over {0!s} longitudes'.format(nave))
    x11.add_cell_method(cm)
    cube_out=x11
    
    return cube_out
        
#==========================================================================

def standardise_time_coord_units(cube,timename='time',tunits=False,basetime=' since 1900-01-01 00:00:0.0',verbose=True):
    """Convert time coordinate of cube to standard time units.

    Takes account of calendar of time units.

    Standard time units are '[days,hours,seconds] since 1900-01-01 00:00:0.0
    (days,hours, or seconds depending on time units of input cube), or can be
    overwritten with tunits argument and/or basetime argument.
    """
    tcoord1=cube.coord(timename)
    tunits1=tcoord1.units
    #
    if tunits:
        if tunits not in ['days','hours','seconds']:
            raise UserWarning('Invalid tunits.')
        xx=tunits
    else:
        xx=str(tunits1).split()[0]
    tunits2str=xx+basetime
    tunits2=cf_units.Unit(tunits2str,calendar=tunits1.calendar)
    #
    timevals2=tunits1.convert(tcoord1.points,tunits2)
    tcoord2=iris.coords.DimCoord(timevals2,standard_name='time',var_name='time',units=tunits2)
    #
    cube.remove_coord(timename)
    cube.add_dim_coord(tcoord2,0)
    #
    if verbose:
        print('Changing time coordinate to standard base time.')
        print('tunits1: {0!s}'.format(tunits1.__repr__()))
        print('tunits2: {0!s}'.format(tunits2.__repr__()))
        print('tcoord1 [0],[-1]: {0!s},{1!s}'.format(tcoord1.points[0],tcoord1.points[-1]))
        print('tcoord2 [0],[-1]: {0!s},{1!s}'.format(tcoord2.points[0],tcoord2.points[-1]))
    #
    return cube

#==========================================================================

def check_coord_values_match(cube1,cube2):
    """Check the coordinate values in two cubes match.

    Inputs: <cube1>,<cube2>: Two iris cubes.

    Checks that both cubes have the same coordinate names, and that
    all the coordinate values are identical.

    If they do, return True.
    If they do not, raise an exception.
    """

    coordnames=[xx.name() for xx in cube1.coords()]
    coordnames2=[xx.name() for xx in cube2.coords()]
    if coordnames!=coordnames2:
        raise UserWarning('Coordinate names do not match between cubes.')
    for coordnamec in coordnames:
        print('Checking all values in {0!s} coordinate match.'.format(coordnamec))
        if np.all(cube1.coord(coordnamec).points==cube2.coord(coordnamec).points):
            print('Ok. All values match.')
        else:
            raise UserWarning('Values do not match.')
    return True

#==========================================================================

def archive_file(aa,filein):
    """Archive file.
    
    Inputs:

    <aa> : object with string attributes basedir and basedir_archive

    <filein> : string of filename (full pathname) to be archived. filein string must begin with aa.basedir

    Copy (archive) filein to same file structure, but under aa.basedir_archive
    """
    print('Archiving.')
    nchar=len(aa.basedir)
    if filein[:nchar]!=aa.basedir:
        raise UserWarning('filein must begin with basedir.')
    file_archive=filein.replace(aa.basedir,aa.basedir_archive)
    print('filein: {0!s}'.format(filein))
    print('file_archive: {0!s}'.format(file_archive))
    shutil.copyfile(filein,file_archive)

#==========================================================================

def f_wind_from_direction(uwnd,vwnd):
    """Calculate wind_from_direction.

    Inputs: uwnd, vwnd : iris cubes of same dimension

    Output: wnddir : iris cube of same dimension as uwnd

    wnddir is wind direction as a bearing in degrees, giving the direction the wind is
    coming FROM, as in the standard meteorological convention.
    
    Examples
    uwnd vwnd wnddir   name
      0    3     180   southerly
      3    3     225   southwesterly
      3    0     270   westerly
      3   -3     315   northwesterly
      0   -3       0   northerly
     -3   -3      45   northeasterly
     -3    0      90   easterly
     -3    3     135   southeasterly
    """
    x1=(360/math.tau)*np.arctan2(-vwnd.data,-uwnd.data)
    x2=90-x1
    x3=np.where(x2<0,x2+360,x2)
    wnddir=create_cube(x3,uwnd,new_var_name='wnddir')
    return wnddir

#==========================================================================

def f_south2north(cube,verbose=False):
    """Determine whether latitude axis of iris cube runs south2north or not.

    Input: cube : iris cube

    Output: south2north : Boolean variable. True if latitude axis of
    cube runs south to north, False if it runs north to south.

    Return south2north.
    """
    lat_coord=cube.coord('latitude')
    lat_beg=lat_coord.points[0]
    lat_end=lat_coord.points[-1]
    if lat_beg<lat_end:
        south2north=True
    else:
        south2north=False
    if verbose:
        print('f_south2north: lat_beg, lat_end, south2north: {0!s}, {1!s}, {2!s} '.format(lat_beg,lat_end,south2north))
    return south2north

#==========================================================================

def add_month(dt,verbose=False):
    """Add month to datetime and return new datetime.datetime object."""
    yearc=dt.year
    monthc=dt.month
    dayc=dt.day
    hourc=dt.hour
    minutec=dt.minute
    secondc=dt.second
    monthnew=monthc+1
    if monthnew==13:
        monthnew=1
        yearnew=yearc+1
    else:
        yearnew=yearc
    dtnew=datetime.datetime(yearnew,monthnew,dayc,hourc,minutec,secondc)
    if verbose:
        print('# add_month')
        print(dt)
        print(dtnew)
    return dtnew

#==========================================================================

def f_diurnal_cycle_cube(cube_in,source,verbose=True):
    """Calculate mean diurnal cycle from iris cube.

    Input: cube_in : iris cube

    cube_in must have a 'time' coordinate that is dimension 0, and no
    missing data or missing times.

    If cube_in does not have an integer number of days, calculation
    will be carried out on data from whole days only, discarding last
    partial day of data in cube_in.

    Output: cube_dc : iris cube of mean diurnal cycle calculcated from cube_in.

    Return cube_dc.

    """
    # Check cube_in has time coordinate at dimension 0.
    if cube_in.dim_coords[0].name()!='time':
        raise UserWarning('First dimension of cube needs to be time.')
    # Extract integer number of days from cube_in
    npd=find_npd(source)
    tcoord=cube_in.coord('time')
    tunits=tcoord.units
    ntime=len(tcoord.points)
    print('ntime,npd: {0!s},{1!s}'.format(ntime,npd))
    nfulldays=divmod(ntime,npd)[0]
    ntimeextract=nfulldays*npd
    if divmod(ntime,npd)[1]==0:
        print('cube_in has an integer number of days.')
        x1=cube_in
    else:
        print('cube_in does not have an integer number of days.')
        print('nfulldays,ntimeextract: {0!s},{1!s}'.format(nfulldays,ntimeextract))
        x1=cube_in[:ntimeextract]
    # Calculate mean for each time of day
    x2=iris.cube.CubeList([])
    tunitsdc=cf_units.Unit('minutes since 1000-01-01')
    for ii in range(npd):
        print('ii: {0!s}'.format(ii))
        # Extract data for current time of day
        x3=x1[ii:ntimeextract:npd]
        tcoordii=x3.coord('time')
        print(tcoordii)
        x4=x3.collapsed('time',iris.analysis.MEAN)
        # Set time axis for diurnal cycle contribution to 1 Jan 1000
        time1=tunits.num2date(tcoordii.points[0])
        timec=datetime.datetime(1000,1,1,time1.hour,time1.minute,time1.second)
        timecval=tunitsdc.date2num(timec)
        tcoordc=iris.coords.DimCoord([timecval],standard_name='time',units=tunitsdc)
        x5=create_cube(x4.data.reshape((1,)+x4.shape),x1,new_axis=tcoordc)
        x2.append(x5)
        #hhhhhhh
    cube_dc=x2.concatenate_cube()
    # Add cell method to describe calculation of diurnal cycle
    cm=iris.coords.CellMethod('point','time',comments='mean diurnal cycle')
    cube_dc.add_cell_method(cm)
    #
    return cube_dc

#==========================================================================

def f_subtract_diurnal_cycle_cube(cube_in,cube_dc,source,verbose=True):
    """Subtract mean diurnal cycle from iris cube.

    Inputs:

    cube_in : iris cube. Same restrictions as input to f_diurnal_cycle_cube.

    cube_dc: iris cube, output from f_diurnal_cycle_cube or similar.

    Output: cube_out : iris cube of cube_in with diurnal cycle of
    cube_dc subtracted.
    """
    pdb.set_trace()

#==========================================================================

class ToDoError(UserWarning):

    """An exception that indicates I need to write some more code.

    Use raise ToDoError as a placeholder for code to be written for a
    different situation."""

    pass

#==========================================================================

class Dummy(object):

    """A dummy class.

    Can set any attributes. Use when a class object is needed with certain attributes."""

    def __init__(self):
        """Initialise dummy class object."""
        self.verbose=True

#==========================================================================

class TimeDomain(object):

    """A set of single or paired (start,end) times.

    A set of single times can be used as the basis of (lagged)
    composite analysis.

    A set of paired (start,end) times can be used as the basis of
    calculating time means and other statistics.

    They have a unique string identifier, e.g., 'jul9899'.

    They are stored/archived as ascii files in
    /gpfs/home/e058/home/data/tdomain/

    Once created, the ascii time domain files should not be modified
    or deleted.  A modified version of the time domain should be given
    a new id, and a new ascii file created for it.

    For historical reasons (from using the cdat module) the format of
    the ascii file is that from a print cdtime.comptime object,
    e.g. for 'jul9899':

    # Optional descriptive comment line(s) beginning with '#'
    # July 1998, 1999
    1998-07-01 00:00:0.0, 1998-07-31 23:59:0.0
    1999-07-01 00:00:0.0, 1999-07-31 23:59:0.0

    This format must be adhered to.

    The TimeDomain class has methods to read and write these ascii
    files, and to convert to/from datetime objects, for use with iris
    scripts.

    Attributes:

    self.idx : unique string identifier, e.g., 'jul9899'.  Must be
    supplied when creating instance of TimeDomain class. Must not
    contain an underscore _.

    self.basedir : directory path for timedomain ascii files.  Has
    default value that can be overwritten as a keyword argument.

    self.filename : full path name for timedomain ascii file.  Created
    from self.basedir and self.idx.

    self._format2ascii : format used for conversion to ascii format.
    Do not change.

    self._format2datetime : format used for conversion to datetime
    object.  Do not change.

    self.lines : ascii representation of timedomain.  A nested list of
    strings.

    self.datetimes : datetime representation of timedomain.  A nested
    list of datetime objects.

    self.type : either 'single' or 'event'.

    self.nevents : integer equivalent to length of self.lines, i.e., number
    single times or events.

    self.header : tuple of two strings, each beginning with '#', that
    are descriptors of the time domain

    self.ndaytot : float of total number of days covered by time
    domain if it is of type 'event'.

    self.ndays : sorted list of floats of number of days in each event
    in time domain, if time domain is of type 'event'.

    Note on calendars. The TimeDomain class works with
    datetime.datetime objects which implicitly assume a gregorian
    calendar. When using 360_day calendar data, still use the
    TimeDomain class, and call the set_time_constraint function when
    creating iris time constraints from the TimeDomain
    datetime.datetime objects. This will convert the datetime.datetime
    objects to cftime.Datetime360Day objects needed to extract 360_day
    calendar data.

    However, the actual dates in a time domain might have to change
    according to the calendar. For example, a JJA average time domain
    would have a line like e.g. 1986-06-01 to 1986-08-31 if the
    calendar is gregorian, but 1986-06-01 to 1986-08-30 if the
    calendar is 360_day. Hence, the convention adopted is that time
    domains created for 360_day calendar data have their id ending in
    '-360day', e.g., 'jja86-88-360day'.

    """
    
    def __init__(self,idx,basedir=mypaths.DIR_TDOMAIN_DEFAULT,tseriesdir=mypaths.DIR_TSERIES_DEFAULT,verbose=True):
        self.idx=idx
        self.basedir=basedir
        self.tseriesdir=tseriesdir
        self.verbose=verbose
        self.filename=os.path.join(self.basedir,self.idx+'.txt')
        self._format2ascii="%Y-%m-%d %H:%M"
        self._format2datetime="%Y-%m-%d %H:%M:%S:%f"
        if self.verbose:
            print(self)
        if '_' in self.idx:
            raise UserWarning('Time domain id must not contain an underscore.')

    def __repr__(self):
        return 'TimeDomain({0.idx!r},basedir={0.basedir!r},verbose={0.verbose!r})'.format(self)

    def __str__(self):
        if self.verbose==2:
            ss=h1a+'TimeDomain({0.idx!s},basedir={0.basedir!s},verbose={0.verbose!s}) \n'+\
                'filename={0.filename!s} \n'+h1b
            return ss.format(self)
        else:
            return self.__repr__()

    def read_ascii(self):
        """Read the ascii time domain file.

        Discard any comments (lines beginning with '#')

        Create lines and comments attribute.
        """
        file1=open(self.filename)
        lines=file1.readlines()
        file1.close()
        self.header=[xx for xx in lines if xx[0]=='#']
        self.lines=[xx for xx in lines if xx[0]!='#']
        if self.verbose:
            ss='read_ascii:  {0.idx!s}: Created attributes "comments", "lines". \n'
            if self.verbose==2:
                ss+='   header: {0.header!s} \n'
                ss+='   lines: {0.lines!s} \n'
            print(ss.format(self))

    def write_ascii(self):
        """Write the ascii file for a newly created time domain."""
        if os.path.isfile(self.filename):
            raise UserWarning('Warning: time domain already exists.  Not overwriting.')
        else:
            # Write ascii strings
            file1=open(self.filename,'w')
            if self.header:
                file1.writelines(self.header)
            file1.writelines(self.lines)
            file1.close()
            if self.verbose:
                ss='write_ascii:  {0.idx!s}: Created ascii file. \n'
                print(ss.format(self))

    def ascii2datetime(self):
        """Convert ascii representation of timedomain to datetime representation.
        Create a datetimes attribute.

        """
        datetimes=[]
        for line in self.lines:
            # Separate comma-separated string times
            a=line.strip()
            a=a.split(',')
            datetimes_row=[]
            for t1 in a:
                # Convert decimal seconds, e.g., 0.0 from the ascii format
                # to integer seconds and integer microseconds for the
                # datetime object
                # The seconds will be after the final ":" in the t1 string
                t1a=t1.strip()
                index_decimal_point=0
                found_decimal_point=False
                while not found_decimal_point:
                    index_decimal_point-=1
                    if t1a[index_decimal_point]=='.':
                        found_decimal_point=True
                index_colon=0
                found_colon=False
                while not found_colon:
                    index_colon-=1
                    if t1a[index_colon]==':':
                        found_colon=True
                second='%02d' % int(t1a[index_colon+1:index_decimal_point])
                microsecond='%06d' % (1000000*int(t1a[index_decimal_point+1:]))
                t1b=t1a[:index_colon+1]+second+':'+microsecond
                # Convert string time to a datetime object and append
                xx=datetime.datetime.strptime(t1b,self._format2datetime)
                datetimes_row.append(xx)
            # Append this row of datetime objects to master list
            datetimes.append(datetimes_row)
        self.datetimes=datetimes
        if self.verbose:
            ss='ascii2datetime:  {0.idx!s}: Created attribute "datetimes". \n'
            if self.verbose==2:
                ss+='   datetimes: {0.datetimes!s} \n'
            print(ss.format(self))

    def ascii2partial_date_time(self):
        """Convert ascii representation of timedomain to PartialDateTime.
        
        Create a partial_date_times attribute.

        This method should not be needed.  I was confused over partial
        date times and datetimes.  Just use datetimes.

        """
        raise DeprecationWarning('This method should not be needed.  Use ascii2datetime instead.')
        
        # First create datetimes attribute if it does not exist
        try:
            datetimes=self.datetimes
        except AttributeError:
            self.ascii2datetime()
            datetimes=self.datetimes
        print('datetimes: {0!s}'.format(datetimes))
        # Convert datetime objects to PartialDateTime objects
        partial_date_times=[]
        for datetimes_row in datetimes:
            partial_date_times_row=[]
            for t1 in datetimes_row:
                pdtc=PartialDateTime(year=t1.year,month=t1.month,day=t1.day,hour=t1.hour,minute=t1.minute,second=t1.second,microsecond=t1.microsecond)
                partial_date_times_row.append(pdtc)
            # Append this row of PartialDateTime objects to master list
            partial_date_times.append(partial_date_times_row)
        self.partial_date_times=partial_date_times
        if self.verbose:
            ss='ascii2partial_date_time:  {0.idx!s}: Created attribute "partial_date_times". \n'
            if self.verbose==2:
                ss+='   partial_date_times: {0.partial_date_times!s} \n'
            print(ss.format(self))
        
    def datetime2ascii(self):
        """Convert datetime representation of timedomain to ascii.

        Create a lines attribute.

        """
        # Convert datetime objects to ascii strings
        lines=[]
        for row in self.datetimes:
            lines_row=''
            # Convert datetime object to formatted string
            # Integer seconds and integer microseconds in the datetime
            # object must be converted to decimal seconds for the ascii
            # format
            if isinstance(row,list):
                # 'event' type
                for t1 in row:
                    decimal_second=str(t1.second+float(t1.microsecond)/1e6)
                    xx=datetime.datetime.strftime(t1,self._format2ascii)+':'+decimal_second+', '
                    lines_row+=xx
                # Remove the final ', ' from the last time
                lines_row=lines_row[:-2]
            else:
                # 'single' type
                t1=row
                decimal_second=str(t1.second+float(t1.microsecond)/1e6)
                xx=datetime.datetime.strftime(t1,self._format2ascii)+':'+decimal_second
                lines_row+=xx
            # Add newline character
            lines_row+='\n'
            lines.append(lines_row)
        self.lines=lines
        if self.verbose:
            ss='read_ascii:  {0.idx!s}: Created attribute "lines". \n'
            if self.verbose==2:
                ss+='   lines: {0.lines!s} \n'
            print(ss.format(self))

    def time_domain_type(self):
        """Determine if time domain is 'single' or 'event' type.

        'Single' type is list of single times.

        'Event' type is list of paired (start,end) times.

        Creates type attribute.
        """
        t1=self.datetimes[0]
        if len(t1)==1:
            self.type='single'
        elif len(t1)==2:
            self.type='event'
        else:
            raise UserWarning('Error: There must be either 1 or 2 times in each row.')
        if self.verbose:
            ss='time_domain_type:  {0.idx!s}: Created attribute "type". \n'
            if self.verbose==2:
                ss+='   type: {0.type!s} \n'

    def f_nevents(self):
        """Count number of events in time domain.

        Create attribute nevents, which is length of datetimes attribute.
        """
        self.nevents=len(self.datetimes)

    def f_create_time_domain_from_indices(self):
        """Create time domain from algorithm on index time series."""
        if len(self.idx)==10 and self.idx[:3] in ['rmm','bs1','bs2']:
            # Time domains based on RMM or BSISO indices
            # Format is 'PPPXXXYYYZ' where
            #   'PPP' is name of index: 'rmm','bs1','bs2'
            #   'XXX' is a counter e.g. '001'
            #   'YYY' is a season e.g. 'djf', 'jja', 'n2a', 'm2o'
            #   'Z' is the category ('1' to '8')
            # e.g., rmm001djf1
            indexname=self.idx[0:3]
            counter=self.idx[3:6]
            season=self.idx[6:9]
            category=int(self.idx[9])
            print('indexname: {0!s}'.format(indexname))
            print('counter: {0!s}'.format(counter))
            print('season: {0!s}'.format(season))
            print('category: {0!s}'.format(category))
            # Counters:
            header1='# '+indexname+': '
            if counter=='001':
                header1+='Index amplitude >=1, time range 1 Jan 1979 to 31 Dec 2015 \n'
                time1=datetime.datetime(1979,1,1)
                time2=datetime.datetime(2015,12,31)
                time_constraint=iris.Constraint(time=lambda cell: time1<=cell<=time2)
                amp_threshold=1
            elif counter=='002':
                header1+='Index amplitude >=0.75, time range 1 Jan 1979 to 31 Dec 2015 \n'
                time1=datetime.datetime(1979,1,1)
                time2=datetime.datetime(2015,12,31)
                time_constraint=iris.Constraint(time=lambda cell: time1<=cell<=time2)
                amp_threshold=0.75
            elif counter=='003':
                # Same conditions as Fig. 9 in Lee et al. (2013)
                header1+='Index amplitude >=1.5, time range 1 Jan 1981 to 31 Dec 2010 \n'
                time1=datetime.datetime(1981,1,1)
                time2=datetime.datetime(2010,12,31)
                time_constraint=iris.Constraint(time=lambda cell: time1<=cell<=time2)
                amp_threshold=1.5
            elif counter=='004':
                # As '001', but only from 1980-2014 because of loss
                # of ends of data set from time filtering
                header1+='Index amplitude >=1, time range 1 Jan 1980 to 31 Dec 2014 \n'
                time1=datetime.datetime(1980,1,1)
                time2=datetime.datetime(2014,12,31)
                time_constraint=iris.Constraint(time=lambda cell: time1<=cell<=time2)
                amp_threshold=1
            elif counter=='005':
                # As '001', but for 1998-2018
                header1+='Index amplitude >=1, time range 1 Jan 1998 to 31 Dec 2018 \n'
                time1=datetime.datetime(1998,1,1)
                time2=datetime.datetime(2018,12,31)
                time_constraint=iris.Constraint(time=lambda cell: time1<=cell<=time2)
                amp_threshold=1
            elif counter=='006':
                # As '005', but for 1999-2017
                header1+='Index amplitude >=1, time range 1 Jan 1998 to 31 Dec 2018 \n'
                time1=datetime.datetime(1999,1,1)
                time2=datetime.datetime(2017,12,31)
                time_constraint=iris.Constraint(time=lambda cell: time1<=cell<=time2)
                amp_threshold=1
            else:
                raise UserWarning('Counter is not valid.')
            # Seasons
            if season=='djf':
                valid_months=[12,1,2]
            elif season=='n2a':
                valid_months=[11,12,1,2,3,4]
            elif season=='m2o':
                valid_months=[5,6,7,8,9,10]
            elif season=='all':
                valid_months=[1,2,3,4,5,6,7,8,9,10,11,12]
            else:
                raise UserWarning('season is not valid.')
            # Category
            if category not in list(range(1,8+1)):
                raise UserWarning('category is not valid.')
            header2='# Season='+season+', category='+str(category)+' \n'
            # Set header attribute for time domain
            self.header=(header1,header2)
            print(header1)
            print(header2)
            # Read input time series of index amplitude and category
            if indexname=='rmm':
                source='rmm_nl_d'
                f1=os.path.join(self.tseriesdir,source,'std','rmm_amp_-999.nc')
                f2=os.path.join(self.tseriesdir,source,'std','rmm_cat_-999.nc')
                x1=iris.load_cube(f1,'RMM_amplitude')
                x2=iris.load_cube(f2,'RMM_category')
            elif indexname=='bs1':
                source='bsiso_nl_d'
                f1=os.path.join(self.tseriesdir,source,'std','bsiso1_amp_-999.nc')
                f2=os.path.join(self.tseriesdir,source,'std','bsiso1_cat_-999.nc')
                x1=iris.load_cube(f1,'BSISO_1_amplitude')
                x2=iris.load_cube(f2,'BSISO_1_category')
            else:
                raise UserWarning('Invalid indexname.')
            # Extract data for valid time range
            x1a=x1.extract(time_constraint)
            x2a=x2.extract(time_constraint)
            time_coord=x1a.coord('time')
            time_units=time_coord.units
            # Flags indicate whether relevant time is within an event (True)
            # or not (False)
            flag_previous_time=False
            flag_current_time=False
            datetimes_list=[]
            # Create list of datetime.datetime pairs of (start,end) dates
            #   for events
            # Loop over time
            for timevalc in time_coord.points:
                timecompc=time_units.num2date(timevalc)
                time_constraint2=set_time_constraint(timecompc,False)
                # Extract index amplitude and category at current time
                ampc=float(x1a.extract(time_constraint2).data)
                catc=float(x2a.extract(time_constraint2).data)
                # Set flag for current time
                if ampc>=amp_threshold and catc==category and timecompc.month in valid_months:
                    flag_current_time=True
                else:
                    flag_current_time=False
                print(timevalc,timecompc,ampc,catc,flag_previous_time,flag_current_time)
                # Use flags to determine start and end points of events
                if not flag_previous_time and flag_current_time:
                    print('Start of event')
                    start_event=timecompc
                if flag_previous_time and not flag_current_time:
                    print('End of event')
                    end_event=timecomp_previous_time
                    datetimes_list.append([start_event,end_event])
                # Current time becomes previous time in next iteration
                flag_previous_time=flag_current_time
                timecomp_previous_time=timecompc
            # Set datetimes attribute and write time domain
            self.datetimes=datetimes_list
            self.datetime2ascii()
            self.write_ascii()
        elif self.idx[0]=='M' and len(self.idx)==5:
            # 'Mxxxx' time domains
            # Create a time domain from a pre-existing time series according to the threshold information
            # Time domain id has form Mxxxx where xxxx is an integer (starting at 0001) and the information
            # that defines Mxxxx is in info.py
            dictc=info.Mtdomain[self.idx[1:]]
            self.source=dictc['source']
            tseriesfile=dictc['tseriesfile']
            method=dictc['method']
            print('source: {0.source!s}'.format(self))
            print('tseriesfile: {0!s}'.format(tseriesfile))
            print('method: {0!s}'.format(method))
            source_info(self)
            # Read in time series (probably from multiple files)
            f1=os.path.join(self.tseriesdir,self.source,'processed',tseriesfile+'_'+self.wildcard+'.nc')
            #f1=os.path.join(self.tseriesdir,self.source,'processed',tseriesfile+'_'+'199804'+'.nc')
            print('f1: {0!s}'.format(f1))
            x1=iris.load(f1)
            x1=x1.concatenate_cube()
            print('Successfully read time series and concatenated.')
            # Find threshold value
            if method=='percentile_above':
                percentile=dictc['percentile']
                valc=np.percentile(x1.data,percentile)
                print('{0!s}th percentile is: {1!s}'.format(percentile,valc))
                unitsc=str(x1.units)
                time_coord=x1.coord('time')
                time_units=time_coord.units
                t1comp=time_units.num2date(time_coord.points[0])
                t2comp=time_units.num2date(time_coord.points[-1])
                header1='# Events exceeding '+str(percentile)+'th percentile ('+str(valc)+' '+str(unitsc)+')\n'
                header2='# Taken from '+self.source+' '+tseriesfile+': '+str(t1comp)+' to '+str(t2comp)+'\n'
                self.header=(header1,header2)
                print(header1)
                print(header2)
                x1=np.where(np.greater_equal(x1.data,valc),time_coord.points,0)
                x2=x1[x1!=0]
                self.datetimes=[time_units.num2date(xx) for xx in x2]
                self.datetime2ascii()
                self.write_ascii()
            else:
                raise ToDoError('Code up for other methods.')


    def f_ndays(self):
        """Calculate number of days covered by time domain.

        Create ndays and ndaytot attributes.
        """
        if self.type=='event':
            ndaytot=0
            ndays=[]
            for (time_beg,time_end) in self.datetimes:
                timedeltac=time_end-time_beg
                ndayc=timedeltac.total_seconds()/86400+1
                ndays.append(ndayc)
                ndaytot+=ndayc
            ndays.sort()
            self.ndays=ndays
            self.ndaytot=ndaytot
        else:
            raise UserWarning('Time domain must be of type "event" to calculate number of days.')

    def f_histogram_months(self):
        """Create histogram of occurence of months in (timebeg) of time domain.

        Print statistics of occurence during standard seasons.

        Create histogram_months attribute.
        """
        hist=np.zeros((12,))
        for xx in self.datetimes:
            monthc=xx[0].month
            indexc=monthc-1
            hist[indexc]=hist[indexc]+1
        print('hist: {0!s}'.format(hist))
        ntot=hist.sum()
        print('ntot: {0!s}'.format(ntot))
        ndjf=hist[11]+hist[0]+hist[1]
        fracdjf=ndjf/ntot
        print('ndjf,fracdjf: {0!s}, {1!s}'.format(ndjf,fracdjf))
        nmam=hist[2]+hist[3]+hist[4]
        fracmam=nmam/ntot
        print('nmam,fracmam: {0!s}, {1!s}'.format(nmam,fracmam))
        njja=hist[5]+hist[6]+hist[7]
        fracjja=njja/ntot
        print('njja,fracjja: {0!s}, {1!s}'.format(njja,fracjja))
        nson=hist[8]+hist[9]+hist[10]
        fracson=nson/ntot
        print('nson,fracson: {0!s}, {1!s}'.format(nson,fracson))
        nndjfma=hist[10]+hist[11]+hist[0]+hist[1]+hist[2]+hist[3]
        fracndjfma=nndjfma/ntot
        print('nndjfma,fracndjfma: {0!s}, {1!s}'.format(nndjfma,fracndjfma))
        nmjjaso=ntot-nndjfma
        fracmjjaso=1-fracndjfma
        print('nmjjaso,fracmjjaso: {0!s}, {1!s}'.format(nmjjaso,fracmjjaso))
        self.histogram_months=hist

    def info(self):
        """Calculate and print information on time domain."""
        self.read_ascii()
        self.ascii2datetime()
        self.f_nevents()
        print('nevents: {0.nevents!s}'.format(self))
        self.time_domain_type()
        print('type: {0.type!s}'.format(self))
        if self.type=='event':
            self.f_ndays()
            print('ndays: {0.ndays!s}'.format(self))
            print('ndaytot: {0.ndaytot!s}'.format(self))

    def f_randomise_times(self,idx_new,max_day_shift,time_first,time_last):
        """Create copy of time domain with randomised times.

        Inputs:

        idx_new : string id of new (randomised) time domain

        max_day_shift : integer (typical value 15)

        time_first, time_last : datetime.datetime objects
        corresponding to the time range of the data this randomised
        time domain will be applied to.
        
        For each random time, the year of the original time is
        replaced by a random year between time_first.year and
        time_last.year, and the day is shifted randomly between
        plus/minus max_day_shift.  This procedure ensures the
        randomised time is at approximately the same time of year
        (i.e., same part of seasonal cycle) as the original time.

        If the time domain type is 'event', the randomisation is
        carried out on the start time of each pair of datetimes, and
        the time difference between the start and end time of each
        pair is preserved.

        Returns tdomain_rand, a new (randomised) TimeDomain object.
        """
        # Loop over datetimes in time domain
        datetimes_rand=[]
        for dtc in self.datetimes:
            # Extract start datetime if type is 'event', or only datetime if type is 'single'
            dt0=dtc[0]
            # If type is 'event' extract end datetime
            if self.type=='event':
                dt1=dtc[1]
                timedelta=dt1-dt0
            # Randomise datetime(s) and check they lie within allowed range
            checked_valid_time_range=False
            while not checked_valid_time_range:
                year_rand=random.randrange(time_first.year,time_last.year+1)
                day_increment_rand=random.randrange(-max_day_shift,max_day_shift+1)
                dt0_rand=datetime.datetime(year_rand,dt0.month,dt0.day,dt0.hour,dt0.minute,dt0.second)+datetime.timedelta(days=day_increment_rand)
                print(year_rand,day_increment_rand,dt0,dt0_rand)
                # If type is 'event' change end datetime consistently to preserve time difference between start and end time of this datetime pair
                if self.type=='event':
                    dt1_rand=dt0_rand+timedelta
                    print(dt1_rand)
                # Check randomised times are within allowable range
                if (self.type=='single' and time_first<=dt0_rand<=time_last) or (self.type=='event' and time_first<=dt0_rand and dt1_rand<=time_last):
                    checked_valid_time_range=True
                    print('Randomised datetimes lie inside allowed range.  Proceed.')
                else:
                    print('Randomised datetimes lie outside allowed range.  Randomise again.')
            if self.type=='event':
                datetimes_rand.append([dt0_rand,dt1_rand])
            elif self.type=='single':
                datetimes_rand.append([dt0_rand])
            else:
                raise UserWarning('Invalid time domain type.')
        # Create new time domain from randomised datetimes
        tdomain_rand=TimeDomain(idx_new)
        tdomain_rand.datetimes=datetimes_rand

        return tdomain_rand

    def f_lagged_time_domain(self,idx_new,timedelta_lag):
        """Create copy of time domain with all time values lagged.

        Inputs:

        idx_new : string id of new (randomised) time domain

        timedelta_lag : datetime.timedelta of the lag

        Returns tdomain_lagged, a new (lagged) TimeDomain object.
        """
        print('idx_new, timedelta_lag: {0!s}, {1!s}'.format(idx_new,timedelta_lag))
        # Loop over datetimes in time domain
        datetimes_lagged=[]
        for dtc in self.datetimes:
            dtc_lagged=[xx+timedelta_lag for xx in dtc]
            datetimes_lagged.append(dtc_lagged)
        # Create new time domain from lagged datetimes
        tdomain_lagged=TimeDomain(idx_new)
        tdomain_lagged.datetimes=datetimes_lagged
        print('Original first event: {0!s}'.format(self.datetimes[0]))
        print('Lagged first event: {0!s}'.format(tdomain_lagged.datetimes[0]))

        return tdomain_lagged

#==========================================================================

class DataConverter(object):

    """Converter to iris-friendly standard format netcdf files.

    Called from preprocess.py.

    Using iris, read in data from different sources.  Convert to
    standard cube format and file format.  Write to new netcdf file.

    Code is flexible to allow for new input data sources.

    Attributes:

    self.basedir : string name of base directory for all data

    self.source : string name of data source, e.g., 'ncepdoe_plev_d',
    'ncepdoe_plev_6h', 'ncepdoe_sfc_6h', 'olrcdr_toa_d' etc.  There
    is a different source for each different combination of source
    data set (e.g., NCEP-DOE reanalysis), type of level (e.g.,
    pressure level), and frequency of input data (e.g., daily).  The
    string should be in the format
    '<data_source>_<level_type>_<frequency>'.  In practice, the source
    attribute is used to set the directory in which the netcdf files
    are stored:
       <self.basedir>/<self.source>/raw/ for original data in the
          original format downloaded from data source web site.  Netcdf
          files in this directory are the input to DataConverter
       <self.basedir>/<self.source>/std/ for the converted, standardised
          data, i.e., the output of DataConverter
       <self.basedir>/<self.source>/processed/ for any subsequent analysis
          on the data, e.g., time means etc.  DataConverter does not use this
          directory.

    self.data_source : string name of data source, e.g., 'ncepdoe',
    'olrcdr', etc.  Note that a particular data_source implies a
    particular latitude-longitude grid, e.g., 'ncepdoe' is a 73x144
    grid.  This is relevant for file sizes and self.outfile_frequency.

    self.level_type : string name of level type, e.g., 'plev', 'toa',
    'sfc', 'theta', etc.  This is determined by self.source.

    self.frequency : string denoting time frequency of data.  This is
    used in file names.  It is determined by self.source.  One of:
       '30m' 30-min resolution
       'h' hourly
       '3h' 3-hourly
       '6h' 6-hourly
       'd' daily
       '5d' pentad (5-day)
       '7d' weekly (7-day)

    self.outfile_frequency : string denoting the time coverage of the
    output netcdf file(s).  It is determined by self.source.  One of:
       'year' for separate files for each year, e.g., 1979, 1980
       'month' for separate files for each calendar month, e.g., 197901, 197902
    
    self.var_name : string name of variable, e.g., 'uwnd'.

    self.raw_name : string name of variable in the raw input file.

    self.level : integer single level, e.g., 200 for 200 hPa.  Level naming
    convention:
       1000, 850, 200, etc. for pressure level data (hPa)
       350, 360, etc. for theta level data
       0 for top of atmosphere, e.g., olr
       1 for surface, e.g., sst
       Note there is no ambiguity between eg a 350 K theta level, and
       a 350 hPa pressure level, because of the self.level_type
       attribute, which is determined by self.source.

    self.cube : iris cube of data.  This is the most important
    attribute of DataConverter. The purpose of DataConverter is to
    convert self.cube to standard form and write it as a netcdf file.

    self.file_mask and self.mask.  Default values are False.  If the
    data is to be masked (e.g., SST data to be masked with a land-sea
    mask), then file_mask is the path name for the mask, and mask is
    the array with the mask.

    """
    
    def __init__(self,**descriptor):
        """Initialise from descriptor dictionary.

        Compulsory keywords: 'verbose','source','var_name','level',
        'basedir'.

        Optional keywords: 'file_mask'.
        """
        self.__dict__.update(descriptor)
        self.descriptor=descriptor
        self.name=var_name2long_name[self.var_name]
        source_info(self)
        self.file_mask=False # Default value, overwrite if exists
        self.mask=False # Ditto
        if descriptor['file_mask']:
            self.file_mask=os.path.join(self.basedir,self.source,'raw',descriptor['file_mask'])
            x1=iris.load(self.file_mask,callback=clean_callback)
            x2=x1.concatenate_cube()
            self.mask=x2.data # numpy array of mask data
            if len(self.mask.shape)!=3 and self.mask.shape[0]!=1:
                # Mask should be 2-D (eg lat,lon) but with a third time
                # dimension of length 1
                raise UserWarning('Expecting 3-d mask of shape (1,?,?)')
            if self.source=='sstrey_sfc_7d':
                self.mask=1-self.mask # Switch the 1's and 0's
        if self.data_source in ['hadgem2esajhog']:
            self.calendar='360_day'
        else:
            self.calendar='gregorian'
        if self.verbose:
            print(self)
        
    def __repr__(self):
        return 'DataConverter({0.descriptor!r},verbose={0.verbose!r})'.format(self)

    def __str__(self):
        if self.verbose==2:
            ss=h1a+'DataConverter instance \n'+\
                'source: {0.source!s} \n'+\
                'var_name: {0.var_name!s} \n'+\
                'name: {0.name!s} \n'+\
                'level: {0.level!s} \n'
            if self.file_mask:
                ss+='file_mask: {0.file_mask!s} \n'+\
                     'mask.shape: {0.mask.shape!s} \n'
            ss+=h1b
            return ss.format(self)
        else:
            return self.__repr__()

    def read_cube(self):
        """Read cube from raw input file.

        Code is by necessity ad hoc as it caters for many different
        data sources with different input formats.

        """
        # Set time constraint for current time block
        time1,time2=block_times(self,verbose=self.verbose)
        time_constraint=set_time_constraint(time1,time2,calendar=self.calendar,verbose=self.verbose)
        #
        # Set input file name(s)
        if self.source in ['erainterim_plev_6h']:
            # Input files are in format .../YYYY/MM/DD/ggapYYYYMMDDHH00.nc
            self.filein1=os.path.join(self.basedir,self.source,'raw',str(self.year),'??','??','ggap'+str(self.year)+'??????00.nc')
        elif self.source in ['erainterim_sfc_6h']:
            # Input files are in format .../YYYY/MM/DD/ggasYYYYMMDDHH00.nc
            self.filein1=os.path.join(self.basedir,self.source,'raw',str(self.year),'??','??','ggas'+str(self.year)+'??????00.nc')
        elif self.source in ['erainterim_sfc_d','erainterim_plev_d']:
            #self.filein1=os.path.join(self.basedir,self.source,'raw',self.var_name+str(self.level)+'_'+str(self.year)+'_d.nc')
            raise UserWarning('Need to recode.')
        elif self.source in ['era5trp_plev_h','era5plp_plev_h','era5plp_sfc_h','era5bar_sfc_h','era5mcw_sfc_h']:
            self.filein1=os.path.join(self.basedir,self.source,'raw',self.var_name+'_'+str(self.level)+'_'+str(self.year)+str(self.month).zfill(2)+'.nc')
        elif self.source in ['ncepdoe_plev_6h','ncepdoe_plev_d','ncepncar_plev_d']:
            if self.var_name=='ta':
                self.filein1=os.path.join(self.basedir,self.source,'raw','air'+'.'+str(self.year)+'.nc')
            elif self.var_name=='zg':
                self.filein1=os.path.join(self.basedir,self.source,'raw','hgt'+'.'+str(self.year)+'.nc')
            else:
                self.filein1=os.path.join(self.basedir,self.source,'raw',self.var_name+'.'+str(self.year)+'.nc')
        elif self.source in ['ncepdoe_sfc_d']:
            if self.var_name=='psfc':
                self.filein1=os.path.join(self.basedir,self.source,'raw','pres.sfc'+'.'+str(self.year)+'.nc')
            else:
                raise UserWarning('Invalid source.')
        elif self.source in ['ncepdoegg_zlev_d'] and self.level==10:
            if self.var_name in ['uwnd','vwnd']:
                self.filein1=os.path.join(self.basedir,self.source,'raw',self.var_name+'.10m.gauss'+'.'+str(self.year)+'.nc')
            else:
                raise UserWarning('Invalid source.')
        elif self.source in ['ncepncar_sfc_d',]:
            self.filein1=os.path.join(self.basedir,self.source,'raw',self.var_name+'.sig995.'+str(self.year)+'.nc')
        elif self.source in ['olrcdr_toa_d','olrinterp_toa_d']:
            self.filein1=os.path.join(self.basedir,self.source,'raw',self.var_name+'.day.mean.nc')
        elif self.source in ['ostial4nrttrp_sfc_d','ostial4reptrp_sfc_d']:
            self.filein1=os.path.join(self.basedir,self.source,'raw',self.var_name+'_'+str(self.level)+'_'+str(self.year)+'-'+str(self.month).zfill(2)+'-*.nc')
        elif self.source in ['sg579m031oi01_zlev_h','sg534m031oi01_zlev_h','sg532m031oi01_zlev_h','sg620m031oi01_zlev_h','sg613m031oi01_zlev_h',]:
            self.filein1=os.path.join(self.basedir,self.source,'raw','oi_zt_2m3h_SG'+self.source[2:5]+'.nc')
        elif self.source in ['sstrey_sfc_7d',]:
            if 1981<=self.year<=1989:
                self.filein1=os.path.join(self.basedir,self.source,'raw',self.var_name+'.wkmean.1981-1989.nc')
            elif 1990<=self.year:
                self.filein1=os.path.join(self.basedir,self.source,'raw',self.var_name+'.wkmean.1990-present.nc')
            else:
                raise UserWarning('Invalid year')
        elif self.source in ['trmm3b42v7_sfc_3h']:
            # Inconsistent file naming from NASA DISC
            # 1998-1999 and 2011-2016 files end in .7.nc
            # 2000-2010 files end in .7A.nc
            # Use '.7*.nc' to cover both
            self.filein1=os.path.join(self.basedir,self.source,'raw',str(self.year)+str(self.month).zfill(2),'3B42.'+str(self.year)+str(self.month).zfill(2)+'*.7*.nc')
        elif self.source in ['imergplp_sfc_30m','imergmcw_sfc_30m','imergmts_sfc_30m','imergmt2_sfc_30m','imergnpl_sfc_30m','imergnp2_sfc_30m','imergtrm_sfc_30m']:
            self.filein1=os.path.join(self.basedir,self.source,'raw',str(self.year),str(self.month).zfill(2),'3B-HHR.MS.MRG.3IMERG.'+str(self.year)+str(self.month).zfill(2)+'*.nc')
        elif self.source in ['tropflux_sfc_d']:
            if self.var_name=='lhfd':
                self.filein1=os.path.join(self.basedir,self.source,'raw','lhf_tropflux_1d_'+str(self.year)+'.nc')
            elif self.var_name=='shfd':
                self.filein1=os.path.join(self.basedir,self.source,'raw','shf_tropflux_1d_'+str(self.year)+'.nc')
            elif self.var_name=='nhfd':
                self.filein1=os.path.join(self.basedir,self.source,'raw','netflux_tropflux_1d_'+str(self.year)+'.nc')
            elif self.var_name=='swrd':
                self.filein1=os.path.join(self.basedir,self.source,'raw','swr_tropflux_1d_'+str(self.year)+'.nc')
            elif self.var_name=='lwrd':
                self.filein1=os.path.join(self.basedir,self.source,'raw','lwr_tropflux_1d_'+str(self.year)+'.nc')
            elif self.var_name in ['taux','tauy']:
                self.filein1=os.path.join(self.basedir,self.source,'raw',self.var_name+'_tropflux_1d_'+str(self.year)+'.nc')
            else:
                raise UserWarning('var_name not recognised.')
        elif self.source in ['hadgem2esajhog_plev_d']:
            if self.var_name=='zg':
                self.filein1=os.path.join(self.basedir,self.source,'raw',self.var_name+'_day_HadGEM2-ES_piControl_r1i1p1_*.nc')
        else:
            raise UserWarning('Data source not recognised.')
        #
        # Set level constraint (set to False if none)
        if self.data_source in ['erainterim'] and self.level_type=='sfc':
            level_constraint=iris.Constraint(surface=0)
        elif self.data_source in ['erainterim'] and self.level_type=='plev':
            level_constraint=iris.Constraint(p=self.level)
        elif self.data_source in ['ncepdoe','ncepncar'] and self.level_type=='plev':
            level_constraint=iris.Constraint(Level=self.level)
        elif self.data_source in ['ncepdoegg'] and self.level_type=='zlev':
            level_constraint=iris.Constraint(Level=self.level)
        elif self.data_source in ['hadgem2esajhog'] and self.level_type=='plev':
            level_constraint=iris.Constraint(air_pressure=self.level)
        elif self.source in ['ncepdoe_sfc_d','ncepncar_sfc_d','olrcdr_toa_d','olrinterp_toa_d','sg579m031oi01_zlev_h','sg534m031oi01_zlev_h','sg532m031oi01_zlev_h','sg620m031oi01_zlev_h','sg613m031oi01_zlev_h','sstrey_sfc_7d','imergplp_sfc_30m','imergmcw_sfc_30m','imergmts_sfc_30m','imergmt2_sfc_30m','imergnpl_sfc_30m','imergnp2_sfc_30m','imergtrm_sfc_30m','trmm3b42v7_sfc_3h','tropflux_sfc_d','era5trp_plev_h','era5plp_plev_h','era5plp_sfc_h','era5bar_sfc_h','era5mcw_sfc_h','ostial4nrttrp_sfc_d','ostial4reptrp_sfc_d']:
            level_constraint=False
        else:
            raise ToDoError('Set an instruction for level_constraint.')
        #
        # Set raw_name of variable in raw input data to use in loading
        # 2 Mar 2018. Issue with load_cube. Now can only load on
        # long_name attribute if it exists. Ignores raw_name 
        self.raw_name=self.name
        if self.data_source in ['erainterim',]:
            erainterim_raw_names={'div':'D', 'vrt':'VO', 'uwnd':'U', 'vwnd':'V', 'omega':'W', 'phi':'Z', 'ta':'T', 'psfc':'SP' }
            if self.var_name in erainterim_raw_names.keys():
                self.raw_name=erainterim_raw_names[self.var_name]
            else:
                raise UserWarning('var_name not recognised.')
        elif self.data_source in ['era5trp','era5plp'] and self.level_type=='plev':
            if self.var_name=='uwnd':
                self.raw_name='u'
            elif self.var_name=='vwnd':
                self.raw_name='v'
        elif self.data_source in ['era5plp','era5bar','era5mcw'] and self.level_type=='sfc':
            if self.var_name=='uwnd':
                self.raw_name='u10'
            elif self.var_name=='vwnd':
                self.raw_name='v10'
            elif self.var_name=='ta':
                self.raw_name='t2m'
            else:
                raise UserWarning('Variable not recognised.')
        elif self.data_source in ['ncepdoe','ncepncar']:
            if self.var_name in ['uwnd','vwnd','omega']:
                self.raw_name=self.var_name
            elif self.var_name=='ta':
                self.raw_name='air'
            elif self.var_name=='zg':
                self.raw_name='hgt'
            elif self.var_name=='psfc':
                self.raw_name='pres'
        elif self.data_source in ['hadgem2esajhog']:
            if self.var_name in ['zg']:
                self.raw_name=self.var_name
        elif self.data_source in ['olrinterp',]:
            self.raw_name='olr'
        elif self.data_source in ['ostial4nrttrp_sfc_d','ostial4reptrp_sfc_d']:
            self.raw_name='analysed_sst'
        elif self.data_source in ['sg579m031oi01','sg534m031oi01','sg532m031oi01','sg620m031oi01','sg613m031oi01',]:
            if self.var_name=='tsc':
                self.raw_name='cons_temp'
            elif self.var_name=='sa':
                self.raw_name='abs_salin'
        elif self.data_source in ['tropflux',]:
            if self.var_name=='lhfd':
                self.raw_name='lhf'
            elif self.var_name=='shfd':
                self.raw_name='shf'
            elif self.var_name=='nhfd':
                self.raw_name='netflux'
            elif self.var_name=='swrd':
                self.raw_name='swr'
            elif self.var_name=='lwrd':
                self.raw_name='lwr'
            elif self.var_name in ['taux','tauy']:
                self.raw_name=self.var_name
        #
        # Load cube using a constraint on var_name because if there is a
        # long_name attribute in the netcdf file this will take precendence
        # over var_name if just using a standard load_cube call.
        if self.source in ['sstrey_sfc_7d','imergplp_sfc_30m','imergmcw_sfc_30m','imergmts_sfc_30m','imergmt2_sfc_30m','imergnpl_sfc_30m','imergnp2_sfc_30m','imergtrm_sfc_30m','trmm3b42v7_sfc_3h','ncepdoegg_zlev_d','ostial4nrttrp_sfc_d','ostial4reptrp_sfc_d']:
            print('# Constraint does not work with data sources listed')
            if level_constraint:
                xx=iris.load(self.filein1,constraints=level_constraint & time_constraint,callback=clean_callback)
            else:
                xx=iris.load(self.filein1,constraints=time_constraint,callback=clean_callback)
        else:
            print('# Load using var_name')
            var_constraint=iris.Constraint(cube_func=(lambda c: c.var_name==self.raw_name))
            if level_constraint:
                if self.source in ['erainterim_plev_6h','erainterim_sfc_6h']:
                    # time constraint does not work with erainterim, but redundant as file name constrains time
                    xx=iris.load(self.filein1,constraints=var_constraint & level_constraint,callback=clean_callback)
                else:
                    xx=iris.load(self.filein1,constraints=var_constraint & level_constraint & time_constraint,callback=clean_callback)
            else:
                xx=iris.load(self.filein1,constraints=var_constraint & time_constraint,callback=clean_callback)
        #
        # Convert time coordinate to standard for erainterim
        # and rewrite latitude coordinate
        # For vrt 500 (and presumably for other variables) some of the values on the 
        # latitude coordinate change slightly (order 1e-6). There seems to be one
        # set of latitude values before 2012-02-29:1800 and another slightly different
        # set after 2012-03-01:0000. Have not checked this exhaustively, but is true
        # for 2012 data and first day in 2013 at least. Impact of this is cannot use
        # concatenate_cube for data either side of this critical date. Solution is to
        # take the latitude axis of the very first time (1979-01-01:0000), and replace
        # the latitude coordinate of all data with these values.
        if self.source in ['erainterim_plev_6h','erainterim_sfc_6h']:
            for cubec in xx:
                cubec=standardise_time_coord_units(cubec,timename='t',tunits='days')
            xx=xx.concatenate()
            # Extract latitude axis of 1979-01-01:0000 data and overwrite this for all data
            if self.source=='erainterim_plev_6h':
                filei1lat=os.path.join(self.basedir,self.source,'raw','1979','01','01','ggap197901010000.nc')
            elif self.source=='erainterim_sfc_6h':
                filei1lat=os.path.join(self.basedir,self.source,'raw','1979','01','01','ggas197901010000.nc')
            else:
                raise UserWarning('Invalid source.')
            print('filei1lat: {0!s}'.format(filei1lat))
            x88=iris.load_cube(filei1lat,var_constraint & level_constraint,callback=clean_callback)
            latcoord1=x88.coord('latitude')
            # Standardise latitude coordinate
            for cubec in xx:
                cubec.remove_coord('latitude')
                cubec.add_dim_coord(latcoord1,1)
        #
        # Convert time coordinate to standard for ostia
        if self.source in ['ostial4nrttrp_sfc_d','ostial4reptrp_sfc_d']:
            for cubec in xx:
                cubec=standardise_time_coord_units(cubec,tunits='days')
        #
        # Number of cubes in cubelist should be 1 except if source in list below
        if self.source not in ['hadgem2esajhog_plev_d','erainterim_plev_6h','erainterim_sfc_6h','ostial4nrttrp_sfc_d','ostial4reptrp_sfc_d']:
            ncubes=len(xx)
            if ncubes!=1:
                raise UserWarning('Not a single cube. ncubes='+str(ncubes))
        #
        # Concatenate cube list to single cube
        self.cube=xx.concatenate_cube()
        #
        # Hack for new netcdf4 ncepdoe which have physically implausible time bounds
        xx=self.cube.coord('time')
        xx.bounds=None
        #
        # Level conversion if needed
        if self.source in ['hadgem2esajhog_plev_d']:
            self.cube.coord('air_pressure').convert_units('hPa') # Convert Pa to hPa
        #
        # Apply mask if appropriate
        if self.file_mask:
            print('Applying mask')
            if self.cube.dim_coords[0].name()!='time':
                raise UserWarning('First dimension of cube needs to be time.')
            # Broadcast ntime copies of mask (1,?,?) to (ntime,?,?)
            ntime=self.cube.shape[0]
            ngrid=self.mask.shape[1]*self.mask.shape[2]
            x1=self.mask.reshape((1,ngrid))
            ones=np.ones((ntime,1))
            x2=np.dot(ones,x1)
            x3=x2.reshape((ntime,self.mask.shape[1],self.mask.shape[2]))
            # Apply mask
            self.cube.data=np.ma.array(self.cube.data,mask=x3)
        if self.verbose==2:
            ss=h2a+'read_cube. \n'+\
                'time1: {0!s} \n'+\
                'time2: {1!s} \n'+h2b
            print(ss.format(time1,time2))
        
    def format_cube(self):
        """Change cube to standard format.

        There a few standard format changes applied to all data sets,
        followed by changes specific to particular data sets.
        
        """
        #
        # Universal format changes
        self.cube.rename(self.name)
        self.cube.var_name=self.var_name
        self.cube.coord('time').bounds=None
        self.cube.attributes['data_source']=self.data_source
        self.cube.attributes['level_type']=self.level_type
        self.cube.attributes['frequency']=self.frequency
        self.cube.attributes['outfile_frequency']=self.outfile_frequency
        #
        # BoBBLE OI glider data from Ben Webber
        if self.data_source[:2]=='sg' and self.data_source[6:]=='031oi01':
            # Missing data is set to zero.  Change to 1e20 and mask
            missing_value=1e20
            x1=np.where(np.equal(self.cube.data,0),missing_value,self.cube.data)
            self.cube.data=np.ma.masked_equal(x1,missing_value)
            ## Reset time, as original hourly data time axis was in days
            ## since ..., but values were only stored to 2 dp.
            ## As 1 hour = 0.041666667 days, significant round off error
            ## New time axis is in hours since first time
            #tc=self.cube.coord('time')
            #time0_val=tc.points[0]
            #if time0_val!=int(time0_val):
            #    raise UserWarning('Need a integer value to start with.')
            #time0_datetime=tc.units.num2date(time0_val)
            #new_time_units='hours since '+str(time0_datetime)
            #ntime=tc.points.shape[0]
            #print('ntime : {0!s}'.format(ntime))
            #print('time0_datetime : {0!s}'.format(time0_datetime))
            #new_time_vals=np.arange(ntime)
            #new_time_coord=iris.coords.DimCoord(new_time_vals,standard_name='time',units=new_time_units)
            #print('new_time_coord : {0!s}'.format(new_time_coord))
            #self.cube.remove_coord('time')
            #self.cube.add_dim_coord(new_time_coord,0)
            # Dec 2016.  New data set with time axis units of
            # 'seconds since 1970-01-01'.  However, this suffers from
            # round off error.
            # 5 Dec 2016.  Ben redid analysis and set time axis to
            # 'hours since 2016-1-1' with first time value being 4369 and the
            # last being 4800 (ntime=432)
            # 9 Jan 2017.  Ben supplied new OI data for SG532 with corrected
            #   salinity.  However, time axis had reverted back to
            #   'seconds since 1970-1-1' which suffers from round off error.
            # Overwrite
            tc=self.cube.coord('time')
            if tc.points.shape!=(432,):
                raise UserWarning('Time axis not right')
            new_time_vals=conv_float32(np.arange(4369,4800+1))
            new_time_units='hours since 2016-1-1'
            new_time_coord=iris.coords.DimCoord(new_time_vals,standard_name='time',units=new_time_units)
            print('new_time_coord : {0!s}'.format(new_time_coord))
            self.cube.remove_coord('time')
            self.cube.add_dim_coord(new_time_coord,0)
            if self.var_name=='tsc':
                self.cube.units='degC'
            elif self.var_name=='lon':
                self.cube.units='degree_east'
        #
        # Reynolds SST weekly data.
        # Time stamp is at beginning of week.
        # Change so it is at the end of the week, by adding 3 (days).
        # NB This can bump a data point at the end of the year into the next
        # year, eg 1982-12-29 is changed to 1984-01-01.  This should not
        # matter as the next step with this data is to linearly interpolate
        # to daily data, which uses all data, not the individual yearly files.
        if self.source=='sstrey_sfc_7d':
            tcoord=self.cube.coord('time')
            time_units=tcoord.units
            if 'day' not in time_units.name:
                raise UserWarning('Expecting time units in days.')
            time_val=tcoord.points+3
            tcoord2=iris.coords.DimCoord(time_val,standard_name='time',units=time_units)
            if self.verbose==2:
                print('tcoord: {0!s}'.format(tcoord))
                print('tcoord2: {0!s}'.format(tcoord2))
            self.cube.remove_coord('time')
            self.cube.add_dim_coord(tcoord2,0)
            print('Added 3 days to time coord so time is now in centre of 7-day mean.')
        #
        # Tropflux daily data.
        # Time stamp is at 12 UTC. Change to 00 UTC by subtracting 0.5 (days).
        # Also, longitude runs from 30.5 to 279.5 (19.5E), missing 10 deg long over Africa
        # Set it to run from 0.5 to 359.5 with no missing longitude
        if self.source=='tropflux_sfc_d':
            # Time stamp
            change_time_stamp_from_12_to_00(self,verbose=self.verbose)
            # Longitude
            # 0.5-19.5E
            lon_constraint1=iris.Constraint(longitude = lambda cell: 360.5 <= cell <= 379.5)
            x1=self.cube.extract(lon_constraint1)
            # create cube of missing values at missing longitudes 20.5-29.5E
            shape=x1.shape[:-1]+(10,)
            x2=1e20*conv_float32(np.ones(shape))
            x2=np.ma.array(x2,mask=x2)
            # 30.5-359.5E
            lon_constraint3=iris.Constraint(longitude = lambda cell:  30.5 <= cell <= 359.5)
            x3=self.cube.extract(lon_constraint3)
            # Concatenate the three longitude bands
            x4=np.ma.concatenate((x1.data,x2,x3.data),2)
            # Create new longitude axis
            lon_val=conv_float32(np.arange(0.5,359.5+1e-6,1.0))
            lon_units=x1.coord('longitude').units
            loncoord=iris.coords.DimCoord(lon_val,standard_name='longitude',units=lon_units,circular=True)
            # Create iris cube
            kdim=0
            dim_coords=[]
            for xx in x1.dim_coords[:-1]:
                dim_coords.append([xx,kdim])
                kdim+=1
            dim_coords.append([loncoord,kdim])
            x5=iris.cube.Cube(x4,units=x1.units,attributes=x1.attributes,cell_methods=x1.cell_methods,dim_coords_and_dims=dim_coords)
            x5.rename(x1.name())
            x5.var_name=x1.var_name
            self.cube=x5
        #
        # trmm3b42v7_sfc_3h first data is at 1998-01-01: 03 UTC, not 00 UTC
        # Make a copy of 03 UTC data and set to 00 UTC.
        if self.source=='trmm3b42v7_sfc_3h':
            tcoord=self.cube.coord('time')
            time_units=tcoord.units
            t1=time_units.num2date(tcoord.points[0])
            if t1==datetime.datetime(1998,1,1,3):
                print('trmm3b42v7_sfc_3h. Copying 1998-01-01: 03 data to 1998-01-01: 00.')
                time_constraint=iris.Constraint(time=t1)
                x1=self.cube.extract(time_constraint)
                x2=iris.util.new_axis(x1,'time') # promote time to dim coord
                x2.remove_coord('time')
                t_new=time_units.date2num(datetime.datetime(1998,1,1,0))
                tcoord_new=iris.coords.DimCoord([t_new,],standard_name='time',var_name='time',units=time_units)
                x2.add_dim_coord(tcoord_new,0)
                x3=iris.cube.CubeList([x2,self.cube])
                x4=x3.concatenate_cube()
                self.cube=x4
        #
        # erainterim_sfc_[d,6h] remove spurious vertical coordinate 'surface'
        # which is set to value 0
        if self.data_source=='erainterim' and self.level_type=='sfc':
            self.cube.remove_coord('surface')
        #
        # hadgem2es daily data.
        # Time stamp is at 12 UTC. Change to 00 UTC by subtracting 0.5 (days).
        if self.source in ['hadgem2esajhog_plev_d']:
            change_time_stamp_from_12_to_00(self,verbose=self.verbose)
        #
        # MetUM GOML data.
        # Each year of data has a different time base unit.
        # All the data for 2008 is 'days since 2008-01-01'
        # All the data for 2009 is 'days since 2009-01-01', etc.
        # This means cubes cannot be concatenated.
        # Change the time coordinate to a standard base unit.
        if self.source in ['metumgomlu-bd818_sfc_d']:
            self.cube=standardise_time_coord_units(self.cube,verbose=self.verbose)
        #
        # OSTIA daily data
        # Time stamp is at 12 UTC. Change to 00 UTC by subtracting 0.5 (days).
        if self.source in ['ostial4nrttrp_sfc_d','ostial4reptrp_sfc_d']:
            # Time stamp
            change_time_stamp_from_12_to_00(self,verbose=self.verbose)
        #
        #
        # Final step. Convert cube data to 'single' precision for saving
        self.cube=create_cube(conv_float32(self.cube.data),self.cube)

    def write_cube(self):
        """Write standardised iris cube to netcdf file."""
        # Set output file name
        if self.outfile_frequency=='year':
            self.fileout1=os.path.join(self.basedir,self.source,'std/',self.var_name+'_'+str(self.level)+'_'+str(self.year)+'.nc')
        elif self.outfile_frequency=='month':
            self.fileout1=os.path.join(self.basedir,self.source,'std/',self.var_name+'_'+str(self.level)+'_'+str(self.year)+str(self.month).zfill(2)+'.nc')
        else:
            raise ToDoError("Need to write code for this outfile_frequency.")
        # Write cube
        iris.save(self.cube,self.fileout1)
        if self.archive:
            archive_file(self,self.fileout1)
        if self.verbose==2:
            print('write_cube: {0.fileout1!s}'.format(self))

#==========================================================================

class TimeDomStats(object):

    """Time mean and other statistics of data over non-contiguous time domain.

    Called from several scripts, including mean.py.

    All statistics calculations will be done by this class, so add
    further functions as needed.

    Selected attributes:

    self.tdomainid : string id of the time domain.

    self.tdomain : instance of TimeDomain class for the time domain.

    self.filein1 : path name for input files for data to calculate
    statistics from.  Will likely contain wild card characters.

    self.fileout1 : path name for output file of calculated statistic.
    
    self.cube_event_means : iris cube list of individual event means,
    i.e., the time mean over each pair of (start,end) times in the
    time domain.

    self.cube_event_ntimes : list of integers with number of times
    that went into each even mean, e.g., list of numbers of days in
    each event is data is daily.

    self.time_mean : iris cube of time mean calculated over the whole
    time domain.

    # Following attributes are used in Monte Carlo simulation to
    # calculate null distribution of mean.

    self.nmc : integer number of Monte Carlo simulations to be
    performed for calculation of null distribution.

    self.percentiles_null : list of the percentiles of the null
    distribution to be saved.

    self.max_day_shift : integer maximum day shift to be used in
    randomisation of time domain times for use in calculation of null
    distribution.

    self.time_first : datetime.datetime object that is first (i.e.,
    lowest) allowable date for the randomised times in the Monte Carlo
    simulation.

    self.time_last : as self.time_first but the last (i.e., highest)
    allowable date.

    """

    def __init__(self,lazy_load=True,**descriptor):
        """Initialise from descriptor dictionary.

        Compulsory keywords: 'verbose','source','var_name','level',
        'basedir','tdomainid','filepre','data_from_anncycle'.

        Optional keywords: 'nmc','percentiles_null','max_day_shift',
        'time_first','time_last'.
        """
        self.__dict__.update(descriptor)
        self.descriptor=descriptor
        self.name=var_name2long_name[self.var_name]
        source_info(self)
        # Input data for regular data
        self.filein1=os.path.join(self.basedir,self.source,'std',self.var_name+'_'+str(self.level)+self.filepre+'_'+self.wildcard+'.nc')
        self.data_in=iris.load(self.filein1,self.name)
        if not lazy_load:
            # Force a hard load of entire input data set into memory
            # Potentially speeds up processing
            time_constraint=set_time_constraint(self.time_first,self.time_last,calendar=self.calendar,verbose=self.verbose)
            self.data_in=self.data_in.extract(time_constraint)
        # Input data if calculating statistics from a selection of the annual cycle
        # e.g., to calculate the mean background state over dates in a time domain
        if self.data_from_anncycle:
            file_anncycle=os.path.join(self.basedir,self.source,'processed',self.var_name+'_'+str(self.level)+'_ac_smooth_'+str(self.data_from_anncycle[0])+'_'+str(self.data_from_anncycle[1])+'.nc')
            self.data_in_anncycle=iris.load_cube(file_anncycle,self.name)
            tcoord=self.data_in_anncycle.coord('time')
            tunits=tcoord.units
            self.data_from_anncycle_year=tcoord.units.num2date(tcoord.points[0]).year
        # Time domain
        try:
            dummy1=self.tdomain
            print('tdomain attribute already exists.')
        except AttributeError:
            self.tdomain=TimeDomain(self.tdomainid,verbose=self.verbose)
            self.tdomain.read_ascii()
            self.tdomain.ascii2datetime()
            print('Create tdomain attribute.')
        self.tdomain.time_domain_type()
        self.tdomain.f_nevents()
        # Output files
        self.fileout_mean=os.path.join(self.basedir,self.source,'processed',self.var_name+'_'+str(self.level)+self.filepre+'_'+self.tdomainid+'.nc')
        self.fileout_lagged_mean=os.path.join(self.basedir,self.source,'processed',self.var_name+'_'+str(self.level)+self.filepre+'_'+self.tdomainid+'_lag.nc')
        self.fileout_dc=os.path.join(self.basedir,self.source,'processed',self.var_name+'_'+str(self.level)+self.filepre+'_'+self.tdomainid+'_dc.nc')
        if self.data_from_anncycle:
            # Overwrite output file names (lagged mean, and diurnal cycle files not defined)
            self.fileout_mean=os.path.join(self.basedir,self.source,'processed',self.var_name+'_'+str(self.level)+self.filepre+'_anncycle_'+self.tdomainid+'.nc')
            self.fileout_lagged_mean=''
            self.fileout_dc=''
        if self.verbose:
            print(self)
        
    def __repr__(self):
        return 'TimeDomStats({0.descriptor!r},verbose={0.verbose!r})'.format(self)

    def __str__(self):
        if self.verbose==2:
            ss=h1a+'TimeDomStats instance \n'+\
                'var_name: {0.var_name!s} \n'+\
                'level: {0.level!s} \n'+\
                'source: {0.source!s} \n'+\
                'tdomainid: {0.tdomainid!s} \n'+\
                'data_in: {0.data_in!s} \n'+\
                'filein1: {0.filein1!s} \n'+\
                'fileout_mean: {0.fileout_mean!s} \n'+\
                'fileout_lagged_mean: {0.fileout_lagged_mean!s} \n'+\
                'fileout_dc: {0.fileout_dc!s} \n'+h1b
            return ss.format(self)
        else:
            return 'Statistics of '+self.source+' '+self.var_name+str(self.level)+' over '+self.tdomainid

    def event_means(self):
        """Calculate time mean and ntime for each event in time domain.

        Create cube_event_means and cube_event_ntimes attributes.

        """
        # Loop over events in time domain
        cube_event_means=iris.cube.CubeList([])
        cube_event_ntimes=[]
        for eventc in self.tdomain.datetimes:
            if self.tdomain.type=='event':
                time_beg=eventc[0]
                time_end=eventc[1]
                print('time_beg,time_end: {0!s}, {1!s}'.format(time_beg,time_end))
                time_constraint=set_time_constraint(time_beg,time_end,calendar=self.calendar,verbose=self.verbose)
                if self.data_from_anncycle:
                    raise UserWarning('Event type time domains not suitable for use with annual cycle input data.')
            elif self.tdomain.type=='single':
                timec=eventc[0]
                print('timec: {0!s}'.format(timec))
                if self.data_from_anncycle:
                    # Set time constraint to select month,day of current event from annual cycle
                    if self.calendar=='gregorian' and timec.month==2 and timec.day==29:
                        dayc=28
                    else:
                        dayc=timec.day
                    timec=timec.__class__(self.data_from_anncycle_year,timec.month,dayc,0,0)
                    print('timec: {0!s}'.format(timec))
                time_constraint=set_time_constraint(timec,False,calendar=self.calendar,verbose=self.verbose)
            else:
                raise UserWarning('Invalid tdomain.type')
            if self.data_from_anncycle:
                x2=self.data_in_anncycle.extract(time_constraint)
            else:
                x1=self.data_in.extract(time_constraint)
                x2=concatenate_cube(x1)
            if self.tdomain.type=='event':
                ntime=x2.coord('time').shape[0]
            elif self.tdomain.type=='single':
                ntime=1
            cube_event_ntimes.append(ntime)
            if ntime==1:
                x3=x2
            else:
                x3=x2.collapsed('time',iris.analysis.MEAN)
            cube_event_means.append(x3)
        self.cube_event_means=cube_event_means
        self.cube_event_ntimes=cube_event_ntimes
        self.units=x2.units

    def f_time_mean(self,save=True):
        """Calculate time mean over time domain and optionally save to netcdf.

        Calculate this by a weighted (cube_event_ntimes) mean of the
        cube_event_means.  Hence, each individual time (e.g., day) in the
        original data has equal weighting.

        Create attribute time_mean.

        """
        # Contribution from first event mean
        ntime_total=0
        ntime=self.cube_event_ntimes[0]
        x1=self.cube_event_means[0]*float(ntime)
        ntime_total+=ntime
        # Contribution from remaining events
        if self.tdomain.nevents>1:
            for ievent in range(1,self.tdomain.nevents):
                ntime=self.cube_event_ntimes[ievent]
                x1+=self.cube_event_means[ievent]*float(ntime)
                ntime_total+=ntime
        # Calculate mean
        #time_mean=x1/ntime_total
        # iris bug. Sometimes (eg glider tsc all) the line above fails.
        # Work around.  Make a copy of iris cube, and just access data.
        x2=x1.copy()
        x2.data/=ntime_total
        time_mean=x2
        # Set attributes
        time_mean.rename(self.name)
        time_mean.var_name=self.var_name
        time_mean.units=self.units
        # Add cell method to describe time mean
        cm=iris.coords.CellMethod('point','time',comments='mean over time domain '+self.tdomain.idx)
        time_mean.add_cell_method(cm)
        self.time_mean=time_mean
        if save:
            iris.save(self.time_mean,self.fileout_mean)
            if self.archive:
                archive_file(self,self.fileout_mean)

    def f_time_mean_null_distribution_component(self):
        """Calculate null distribution of the time mean and save.

        Perform a Monte Carlo simulation of size self.nmc to create (a
        partial component of) the null distribution of the time mean.
        For each evaluation in the Monte Carlo simulation, create a
        randomized version of the time domain.
        
        The time mean is then calculated for each of the self.nmc
        randomised time domains, building up a null distribution of
        the mean of length self.nmc.

        Note.  The total number of randomised simulations in the
        overall Monte Carlo simulation should be at least 1000 to
        ensure robustness in the tails of the null distribution.  To
        run these sequentially is too time consuming.  Hence, small
        batches of simulations should be run in parallel.  Use this
        function together called from e.g., mean.py, together with
        run_scripts_bsub.py to do this, looping over a variable called
        NODE_NUMBER.  Each individual simulation runs self.nmc
        realisations, hence the total size of the Monte Carlo
        simulation is NODE_NUMBER*nmc which should be minimum
        approximately 1000.

        Creates attributes: None.
        """
        if self.verbose:
            ss=h1a+'f_time_mean_null_distribution_component \n'+\
                'nmc: {0.nmc!s} \n'+\
                'max_day_shift: {0.max_day_shift!s} \n'+\
                'time_first: {0.time_first!s} \n'+\
                'time_last: {0.time_last!s} \n'+h1b
            print(ss.format(self))
        self.tdomain.time_domain_type()
        zfill_length=len(str(self.nmc))
        # Loop over Monte Carlo simulations and create list of randomised
        # time domains
        tdomains_mc=[]
        for imc in range(self.nmc):
            idxc=self.tdomainid+'_'+str(imc).zfill(zfill_length)
            print('imc,idxc: {0!s}, {1!s}'.format(imc,idxc))
            tdomain_rand=self.tdomain.f_randomise_times(idxc,self.max_day_shift,self.time_first,self.time_last)
            tdomains_mc.append(tdomain_rand)
        self.tdomains_mc=tdomains_mc
        # Loop over Monte Carlo simulations again and calculate
        # time means from each randomised time domain
        x2=iris.cube.CubeList([])
        for imc in range(self.nmc):
            # Create new instance of TimeDomStats object
            print('imc: {0!s}'.format(imc))
            descriptorc=copy.copy(self.descriptor)
            descriptorc['tdomainid']=self.tdomains_mc[imc].idx
            descriptorc['tdomain']=self.tdomains_mc[imc]
            x1=TimeDomStats(**descriptorc)
            # Calculate event means and time mean
            x1.event_means()
            x1.f_time_mean(save=False)
            # Add auxiliary coordinate for simulation_number
            index_mc=self.node_number*self.nmc+imc
            print('index_mc: {0!s}'.format(index_mc))
            mccoord=iris.coords.AuxCoord(float(index_mc),var_name='mc')
            x1.time_mean.add_aux_coord(mccoord)
            # Overwrite cell methods to allow merge later
            x1.time_mean.cell_methods=()
            # Append to cube list
            x2.append(x1.time_mean)
        x3=x2.merge_cube()
        # Save this component of the null distribution
        fileout_mean_mc=os.path.join(self.basedir,self.source,'processed',self.var_name+'_'+str(self.level)+self.filepre+'_'+self.tdomainid+'_mc_'+str(self.node_number).zfill(3)+'.nc')
        print('Saving component of null distribution of mean to {0!s}'.format(fileout_mean_mc))
        iris.save(x3,fileout_mean_mc)
        if self.archive:
            archive_file(self,fileout_mean_nc)

    def f_percentiles_null(self):
        """Calculate percenticles of null distribution of time mean and save.

        Read the null distribution of the time mean previously
        calculated using f_time_mean_null_distribution_component().

        The percentiles self.percentiles_null of this distribution are
        then extracted and saved as self.mean_percentiles_null.

        Creates attributes:

        self.mean_percentiles_null : iris cube containing selected
        percentiles of the null distribution.

        """
        if self.verbose:
            ss=h1a+'f_percentiles_null \n'+\
                'percentiles_null: {0.percentiles_null!s} \n'+h1b
            print(ss.format(self))
        # Read null distribution
        filein_mc=os.path.join(self.basedir,self.source,'processed',self.var_name+'_'+str(self.level)+self.filepre+'_'+self.tdomainid+'_mc_???.nc')
        print('filein_mc: {0!s}'.format(filein_mc))
        x1=iris.load(filein_mc,self.name)
        print('Null distribution - cube list: {0!s}'.format(x1))
        x3=x1.concatenate_cube()
        print('Null distribution - single cube: {0!s}'.format(x3))
        # Calculate percentiles of null distribution of mean
        x4=np.percentile(x3.data,self.percentiles_null,axis=0)
        percoord=iris.coords.DimCoord(np.array(self.percentiles_null),var_name='percentile',units='%')
        x5=create_cube(x4,x3,new_axis=percoord)
        self.mean_percentiles_null=x5
        # Save percentiles
        self.fileout_mean_percentiles_null=os.path.join(self.basedir,self.source,'processed',self.var_name+'_'+str(self.level)+self.filepre+'_'+self.tdomainid+'_percentiles.nc')
        print('Saving percentiles of null distribution of mean to {0.fileout_mean_percentiles_null!s}'.format(self))
        iris.save(self.mean_percentiles_null,self.fileout_mean_percentiles_null)
        if self.archive:
            archive_file(self,self.fileout_mean_percentiles_null)

    def f_lagged_mean(self,method=2,lags=[0,]):
        """Calculate time-lagged means over time domain and save.

        Input:

        method: 1 or 2. Method 1 allows for calculation of a discrete
        set of lags but can be slow with large data sets. Method 2
        calculates all lags between a start lag and an end lag, and is
        order nlag faster. Generally, use method 2.

        lags : if method is 1, lags is a list of integer lags (default
        is [0,]). These are combined with self.timedelta to create the
        lags to calculate the mean over.  E.g., if lags=[-5,0,5] and
        self.timedelta is 1 day, then lagged means will be calculated
        at lags of -5, 0, and 5 days.

        If method is 2, lags is a 2-tuple of (start
        datetime.timedelta, end datetime.timedelta).

        Create attributes:

        self.lags : lags

        self.nlags : length of list of lags

        self.lagged_mean : iris cube of lagged means.  Shape is
        (self.nlags,...)  where ... is shape of input cube.  The 'lag'
        axis is a time axis, relative to an arbitrary reference time
        of 1000-01-01 00:00:0.0

        """
        timelag_units='hours since 1000-01-01 00:00:0.0'
        if method==1:
            # Set lags and nlags attributes
            self.lags=lags
            self.nlags=len(lags)
            self.tdomain.time_domain_type()
            # Loop over lags
            cubelist=iris.cube.CubeList([])
            for lagc in self.lags:
                timedelta_lagc=lagc*self.timedelta
                print('lagc,timedelta_lagc : {0!s}, {1!s}'.format(lagc,timedelta_lagc))
                # Create time domain for this lag
                idxc=self.tdomainid+'-'+str(lagc)
                tdomain_lagc=self.tdomain.f_lagged_time_domain(idxc,timedelta_lagc)
                tdomain_lagc.time_domain_type()
                tdomain_lagc.f_nevents()
                # Calculate time mean for this lagged time domain
                bb=copy.copy(self)
                bb.tdomain=tdomain_lagc
                bb.event_means()
                bb.f_time_mean(save=False)
                xx1=bb.time_mean
                hours_lagc=timedelta_lagc.total_seconds()/3600
                timelag_coord=iris.coords.DimCoord([hours_lagc],standard_name='time',units=timelag_units)
                xx1.add_aux_coord(timelag_coord)
                xx1.cell_methods=None
                cubelist.append(xx1)
            lagged_mean=cubelist.merge_cube()
        elif method==2:
            if self.tdomain.type!='single':
                raise UserWarning("Method 2 only works with time domains of type 'single'.")
            if len(lags)!=2:
                raise UserWarning('For method 2, lags must be a 2-tuple.')
            delta_beg=lags[0]
            delta_end=lags[1]
            kount=0
            # Read lagged data block corresponding to first datetime in time domain
            dtc=self.tdomain.datetimes[0]
            timec=dtc[0]
            timebeg=timec+delta_beg
            timeend=timec+delta_end
            print('timec,timebeg,timeend: {0!s}, {1!s}, {2!s}'.format(timec,timebeg,timeend))
            time_constraint=set_time_constraint(timebeg,timeend,calendar=self.calendar,verbose=self.verbose)
            x1=self.data_in.extract(time_constraint)
            x2=concatenate_cube(x1)
            # Create running sum from first block of data
            x2sum=x2.data
            kount+=1
            # Create lag coordinate using first extracted cube
            tcoord=x2.coord('time')
            self.nlags=len(tcoord.points)
            print('nlags: {0.nlags!s}'.format(self))
            lag_units=cf_units.Unit(timelag_units,calendar=self.calendar)
            lag_first=datetime.datetime(1000,1,1)+delta_beg
            lag_vals=[lag_units.date2num(lag_first+xx*self.timedelta).round(8) for xx in range(self.nlags)]
            lag_coord=iris.coords.DimCoord(lag_vals,standard_name='time',units=timelag_units)
            print('lag_coord: {0!s}'.format(lag_coord))
            # Find time coordinate index for recreation of cube later
            dim_coord_names=[xx.var_name for xx in x2.dim_coords]
            tcoord_index=dim_coord_names.index('time')
            print('tcoord_index: {0!s}'.format(tcoord_index))
            # Loop over remaining times in time domain
            for dtc in self.tdomain.datetimes[1:]:
                # Read in lagged data block corresponding to current datetime in time domain
                timec=dtc[0]
                timebeg=timec+delta_beg
                timeend=timec+delta_end
                print('timec,timebeg,timeend: {0!s}, {1!s}, {2!s}'.format(timec,timebeg,timeend))
                time_constraint=set_time_constraint(timebeg,timeend,calendar=self.calendar,verbose=self.verbose)
                x1=self.data_in.extract(time_constraint)
                x2=concatenate_cube(x1)
                # Calculate contribution to running sum
                x2sum=x2sum+x2.data
                kount+=1
            # Calculate average by dividing by number of events
            print('nevents, kount: {0!s}, {1!s}'.format(self.tdomain.nevents,kount))
            if self.tdomain.nevents!=kount:
                raise UserWarning('Mismatch between nevents and kount.')
            x3=x2sum/kount
            # Create cube of lagged mean
            lagged_mean=create_cube(x3,x2)
            # Replace time axis with lag axis
            lagged_mean.remove_coord('time')
            lagged_mean.add_dim_coord(lag_coord,tcoord_index)
        else:
            raise UserWarning('Invalid method.')
        # Add cell method to describe time mean
        cm=iris.coords.CellMethod('point','time',comments='lagged mean over time domain '+self.tdomain.idx)
        lagged_mean.add_cell_method(cm)
        #
        self.lagged_mean=lagged_mean
        iris.save(self.lagged_mean,self.fileout_lagged_mean)
        if self.archive:
            archive_file(self,self.fileout_lagged_mean)

    def f_diurnal_cycle(self,double=True):
        """Calculate mean diurnal cycle over time domain.

        Time axis for mean diurnal cycle runs for one day (1 Jan in
        year 1, ie 01-01-01).

        if double is True, create a copy of the diurnal cycle for a
        second day, ie 2 Jan year 1.  This is to help in plotting
        later.

        Create attribute:

        self.mean_dc : iris cube of mean diurnal cycle calculated over
        time domain self.tdomain.

        Note: initial attempt at coding this function used
        PartialDateTime to extract data on each day for a given time
        of day (e.g. 0300 UTC) and then average.  This took about
        10-100 times longer to run than method here.
        
        """
        # Extract first day in first event of time domain
        xx=self.tdomain.datetimes[0]
        first_day=xx[0]
        last_day=xx[-1]
        print('first_day,last_day: {0!s}, {1!s}'.format(first_day,last_day))
        time1=datetime.datetime(first_day.year,first_day.month,first_day.day)
        time2=datetime.datetime(first_day.year,first_day.month,first_day.day,23,59,59)
        time_constraint=set_time_constraint(time1,time2,calendar=self.calendar,verbose=self.verbose)
        x1=self.data_in.extract(time_constraint)
        x1=x1.concatenate_cube()
        cube1=x1
        # NB the data attribute of an iris cube is a numpy array if there is
        # no missing data, but a numpy masked array if there is missing data
        # To deal with occasional missing data, convert the data attribute
        # to a numpy masked array, and turn the boolean mask into an array of
        # zeroes and ones.  The kount array then keeps track of the number
        # of non-missing data values contributing to the sum at each spatial
        # grid point
        x1ma=np.ma.array(x1.data)
        x1ma_mask_numerical=np.where(x1ma.mask,0,1)
        x1ma_times_mask=x1ma.data*x1ma_mask_numerical
        data_sum=x1ma_times_mask
        kount=np.zeros(x1.data.shape)
        kount+=x1ma_mask_numerical
        print('time1,time2: {0!s}, {1!s}'.format(time1,time2))
        #print('type(x1.data): {0!s}'.format(type(x1.data)))
        #print(kount.shape,x1ma_mask_numerical.shape)
        #
        # Create time coordinate for final diurnal cycle
        # Extract times of day for data on first day
        tc=x1.coord('time')
        time_units=tc.units
        x2=[tc.units.num2date(xx) for xx in tc.points]
        # Reset these times to year=1,month=1,day=1
        times_datetime=[datetime.datetime(1,1,1,xx.hour,xx.minute,xx.second) for xx in x2]
        times_val=[time_units.date2num(xx) for xx in times_datetime]
        time_coord=iris.coords.DimCoord(times_val,standard_name='time',units=time_units)
        print('times_datetime: {0!s}'.format(times_datetime))
        print('times_val: {0!s}'.format(times_val))
        print('time_coord: {0!s}'.format(time_coord))
        #
        # Loop over (any) remaining days in first event of time domain and add data to data_sum
        timedelta_day=datetime.timedelta(days=1)
        time1+=timedelta_day
        time2+=timedelta_day
        while time1<=last_day:
            raise ToDoError('Recode time constraint to use datetime.datetime')
            time_constraint=set_time_constraint(time1,time2,calendar=self.calendar,verbose=self.verbose)
            x1=self.data_in.extract(time_constraint)
            x1=x1.concatenate_cube()
            x1ma=np.ma.array(x1.data)
            x1ma_mask_numerical=np.where(x1ma.mask,0,1)
            x1ma_times_mask=x1ma.data*x1ma_mask_numerical
            data_sum+=x1ma_times_mask
            kount+=x1ma_mask_numerical
            print('time1,time2: {0!s}, {1!s}'.format(time1,time2))
            #print('type(x1.data): {0!s}'.format(type(x1.data)))
            #print(kount.shape,x1ma_mask_numerical.shape)
            time1+=timedelta_day
            time2+=timedelta_day
        # Loop over (any) remaining events in time domain and add data to data_sum
        for xx in self.tdomain.datetimes[1:]:
            first_day=xx[0]
            last_day=xx[-1]
            print('first_day,last_day: {0!s}, {1!s}'.format(first_day,last_day))
            time1=datetime.datetime(first_day.year,first_day.month,first_day.day)
            time2=datetime.datetime(first_day.year,first_day.month,first_day.day,23,59,59)
            while time1<=last_day:
                raise ToDoError('Recode time constraint to use datetime.datetime')
                time_constraint=set_time_constraint(time1,time2,calendar=self.calendar,verbose=self.verbose)
                x1=self.data_in.extract(time_constraint)
                x1=x1.concatenate_cube()
                x1ma=np.ma.array(x1.data)
                x1ma_mask_numerical=np.where(x1ma.mask,0,1)
                x1ma_times_mask=x1ma.data*x1ma_mask_numerical
                data_sum+=x1ma_times_mask
                kount+=x1ma_mask_numerical
                print('time1,time2: {0!s}, {1!s}'.format(time1,time2))
                #print('type(x1.data): {0!s}'.format(type(x1.data)))
                #print(kount.shape,x1ma_mask_numerical.shape)
                time1+=timedelta_day
                time2+=timedelta_day
        # Divide by kount
        #print('kount: {0!s}'.format(kount))
        print('kount.min,kount.max: {0!s}, {1!s}'.format(kount.min(),kount.max()))
        data_mean=data_sum/kount
        # Create iris cube using metadata from first day of input data
        x10=create_cube(data_mean,cube1)
        x10.remove_coord('time')
        x10.add_dim_coord(time_coord,0)
        # Add cell method to describe diurnal cycle
        cm=iris.coords.CellMethod('point','time',comments='mean diurnal cycle over time domain '+self.tdomain.idx)
        x10.add_cell_method(cm)
        # Create mean_dc attribute
        self.mean_dc=x10
        #
        if double:
            print('Create a double diurnal cycle (ie two identical days)')
            # Create a new time axis for 2 Jan year 1
            times_datetime=[datetime.datetime(1,1,2,xx.hour,xx.minute,xx.second) for xx in x2]
            times_val=[time_units.date2num(xx) for xx in times_datetime]
            time_coord=iris.coords.DimCoord(times_val,standard_name='time',units=time_units)
            print('times_datetime: {0!s}'.format(times_datetime))
            print('times_val: {0!s}'.format(times_val))
            print('time_coord: {0!s}'.format(time_coord))
            # Create copy of cube of diurnal cycle
            x11=self.mean_dc.copy()
            # Apply new time axis
            x11.remove_coord('time')
            x11.add_dim_coord(time_coord,0)
            # Create a cube list of the two diurnal cycles and concatenate
            x12=iris.cube.CubeList([self.mean_dc,x11])
            x13=x12.concatenate_cube()
            # Overwrite mean_dc attribute
            self.mean_dc=x13
            
        # Save diurnal cycle
        iris.save(self.mean_dc,self.fileout_dc)
        if self.archive:
            archive_file(self,self.fileout_dc)
            
#==========================================================================

class TimeFilter(object):
    
    """Time filter using rolling_window method of iris cube.

    Called from filter.py.

    Assumes input data has equally spaced time intervals

    Attributes:

    self.weights : 1-d numpy array of filter weights

    self.nn : integer value of order of filter.

    self.nweights : integer number of weights (length of self.weights
    array.  Must be odd.  Filtered output will be a nweight running
    mean of the input data, using self.weights.

    self.nn : Equal to (nweights-1)/2, e.g., if nweights=61, nn=30.
    Filtered output will be a 61-point weighted mean of the 30 input
    data points before the current time, the current time itself, and
    the 30 intput data points after the current time.

    self.data_out : iris cube of filtered output data.  Length of time
    dimension is typically a convenient block of time, e.g., 1 year
    for daily data.

    self.frequency : string to denote frequency of input (and output)
    data, e.g., 'd' for daily data.

    self.timeout1 : datetime object for start time of self.data_out.
    
    self.timeout2 : datetime object for end time of self.data_out.
    
    self.data_in : iris cube of input data to be filtered.  Length of
    time dimension is length of time dimension of self.data_out +
    2*self.nn (self.nn at the beginning, and self.nn at the end).
    
    self.timein1 : datetime object for start time of self.data_in.
    
    self.timein2 : datetime object for end time of self.data_in.

    self.filein1 : path name for file(s) of input data.
    
    self.fileout1 : path name for file of output (filtered) data.
    
    self.splitblock : The rolling_window method used for filtering
    requires large memory, especially if the number of weights is
    large. If enough memory is not available (job will fail with
    memory error), set splitblock to True. This will split each block
    (of, e.g., one year if outfile_frequency is 'year') up into
    smaller blocks (e.g., months), run the rolling_window on each
    smaller block, then recombine.

    """

    def __init__(self,**descriptor):
        """Initialise from descriptor dictionary.

        Compulsory keywords: 'verbose','source','var_name','level',
        'basedir','filter','file_weights','filepre'
        """
        self.__dict__.update(descriptor)
        self.descriptor=descriptor
        self.name=var_name2long_name[self.var_name]
        source_info(self)
        self.f_weights()
        self.filein1=os.path.join(self.basedir,self.source,'std',self.var_name+'_'+str(self.level)+self.filepre+'_'+self.wildcard+'.nc')
        self.data_in=iris.load(self.filein1,self.name)
        xx=self.source.split('_')
        self.frequency=xx[2]
        if self.frequency=='d':
            self.timedelta=datetime.timedelta(days=self.nn)
        elif self.frequency=='h':
            self.timedelta=datetime.timedelta(hours=self.nn)
        else:
            raise ToDoError('data time interval is not days or hours - need more code!')
        if self.verbose:
            print(self)

    def __repr__(self):
        return 'TimeFilter({0.descriptor!r},verbose={0.verbose!r})'.format(self)

    def __str__(self):
        if self.verbose==2:
            ss=h1a+'TimeFilter instance \n'+\
                'filter: {0.filter!s} \n'+\
                'nn: {0.nn!s} \n'+\
                'nweights: {0.nweights!s} \n'+\
                'weights: {0.weights!s} \n'+\
                'filein1: {0.filein1!s} \n'+\
                'frequency: {0.frequency!s} \n'+\
                'data_in: {0.data_in!s} \n'+h1b
            return(ss.format(self))
        else:
            return 'TimeFilter instance'

    def f_weights(self):
        """Create array of file weights.

        Read in file weights from ASCII text file.  Must be in format:
        Line 0: information on filter (not actually used here)
        Remaining lines: each line has a single float filter weight.

        If nlines is number of lines in the filter file, there are
        nweights=nlines-1 filter weights.

        nweights must be odd.

        nn=(nfilter-1)/2

        Create nweights, nn and weights attributes.

        """
        # Read ASCII filter weights
        f1=open(self.file_weights)
        lines=f1.readlines()
        # Discard first (information) line
        lines2=lines[1:]
        self.nweights=len(lines2)
        # Check nweights is odd
        if divmod(self.nweights,2)[1]!=1:
            raise UserWarning('Error: self.nweights must be odd.')
        self.nn=divmod(self.nweights,2)[0]
        weights=[float(xx) for xx in lines2]
        self.weights=np.array(weights)

    def time_filter(self,subtract=False):
        """Filter using the rolling_window cube method and save data.

        If subtract is False, the filtered data set is just returned.

        If subtract is True, then the filtered data is subtracted from
        the original data, and this subtracted data set is returned.
        This option is generally used to create high-pass filtered
        data.  For example, if the filter weights create an N-running
        mean low pass filter data set, this is then subtracted from
        the original data to create a high-pass filtered data set.

        """
        # Set start and end time of output data
        self.timeout1,self.timeout2=block_times(self,verbose=self.verbose)
        # Calculate start and end time of input data
        self.timein1=self.timeout1-self.timedelta
        self.timein2=self.timeout2+self.timedelta
        # Set output file name
        subtract_string=''
        if subtract:
            subtract_string='minus_'
        if self.outfile_frequency=='year':
            self.fileout1=self.filein1.replace(self.wildcard,subtract_string+self.filter+'_'+str(self.year))
        elif self.outfile_frequency=='month':
            self.fileout1=self.filein1.replace(self.wildcard,subtract_string+self.filter+'_'+str(self.year)+str(self.month).zfill(2))
        else:
            raise UserWarning('outfile_frequency not recognised')
        if self.verbose==2:
            ss=h2a+'timein1: {0.timein1!s} \n'+\
                'timeout1: {0.timeout1!s} \n'+\
                'timeout2: {0.timeout2!s} \n'+\
                'timein2: {0.timein2!s} \n'+\
                'fileout1: {0.fileout1!s} \n'+h2b
            print(ss.format(self))
        # Extract input data
        time_constraint=set_time_constraint(self.timein1,self.timein2,calendar=self.calendar,verbose=self.verbose)
        xx1=self.data_in.extract(time_constraint)
        xx1=xx1.concatenate_cube()
        self.data_current=xx1
        # Apply filter
        if self.splitblock:
            if self.outfile_frequency=='year':
                print('Split current year block into months, filter each month, then concatenate filtered output.')
                # Create dummy class object to call block_times and set attributes for each month
                monthc=Dummy()
                monthc.outfile_frequency='month'
                monthc.calendar=self.calendar
                monthc.year=self.year
                monthc.filtered=iris.cube.CubeList([])
                for imonth in range(1,12+1):
                    print('## imonth: {0!s}'.format(imonth))
                    monthc.month=imonth
                    monthc.timeout1,monthc.timeout2=block_times(monthc,verbose=self.verbose)
                    monthc.timein1=monthc.timeout1-self.timedelta
                    monthc.timein2=monthc.timeout2+self.timedelta
                    if self.verbose==2:
                        ss=h2a+'timein1: {0.timein1!s} \n'+\
                            'timeout1: {0.timeout1!s} \n'+\
                            'timeout2: {0.timeout2!s} \n'+\
                            'timein2: {0.timein2!s} \n'+h2b
                        print(ss.format(monthc))
                    monthc.time_constraint=set_time_constraint(monthc.timein1,monthc.timein2,calendar=monthc.calendar,verbose=monthc.verbose)
                    monthc.data_current=self.data_current.extract(monthc.time_constraint)
                    xx1=monthc.data_current.rolling_window('time',iris.analysis.SUM,self.nweights,weights=self.weights)
                    monthc.filtered.append(xx1)
                xx1=monthc.filtered.concatenate_cube()
            else:
                raise ToDoError('Code up for other outfile frequencies.')
        else:
            xx1=self.data_current.rolling_window('time',iris.analysis.SUM,self.nweights,weights=self.weights)

        # Create a cube from this numpy array
        xx1=create_cube(conv_float32(xx1.data),xx1)
        # Subtract filtered data from original data if required
        if subtract:
            # Re-read input data over same time range as output data
            time_constraint2=set_time_constraint(timeout1,timeout2,calendar=self.calendar,verbose=self.verbose) 
            xx2=self.data_in.extract(time_constraint2)
            xx2=xx2.concatenate_cube()
            # Subtract filtered data
            xx2=xx2.data-xx1.data
            xx1=create_cube(conv_float32(xx2),xx1)
        # Add a cell method to describe the time filter
        cm=iris.coords.CellMethod('mean','time',comments='time filter: '+subtract_string+self.filter)
        xx1.add_cell_method(cm)
        # Set time bounds to None type
        xx1.coord('time').bounds=None
        self.data_out=xx1
        # Save data
        iris.save(self.data_out,self.fileout1)
        if self.archive:
            archive_file(self,self.fileout1)

#==========================================================================

class ModifySource(object):
    """Process data such that it's source attribute has to change.

    Examples:
    1) Time average data to lower time resolution e.g., convert from 
       3-hourly to daily mean.
    2) Time interpolate data to higher time resolution using 
       iris.analysis.interpolate.
    3) Fill in missing data.
    4) Regrid data.

    Selected attributes:

    self.source1 : input source, e.g., trmm3b42v7_sfc_3h 3-hourly data.

    self.source2 : output source, e.g., trmm3b42v7_sfc_d daily data.
    The data_source and level_type parts of self.source1 and
    self.source 2 should be identical.

    self.frequency : the frequency part of self.source2, indicating
    what time resolution the input data is to be converted to, e.g.,
    'd' for daily.

    self.year and self.month : year (and month, depending on value of
    self.outfile_frequency) of input data block.

    self.filein1 : path name for input files from self.source1.
    Probably contains wild card characters.

    self.data_in : iris cube list of all input data

    self.cube_in : input cube of data from self.source1 for the
    current time block.

    self.cube_out : output cube of data for current time block to be
    saved under self.source2.

    self.fileout1 : path name for output file under self.source2.
    """
    
    def __init__(self,**descriptor):
        """Initialise from descriptor dictionary.

        Compulsory keywords: 'verbose','source1','source2','var_name','level',
        'basedir','subdir'.
        """
        self.__dict__.update(descriptor)
        self.descriptor=descriptor
        self.name=var_name2long_name[self.var_name]
        self.source=self.source2
        source_info(self)
        if self.subdir=='std':
            self.filein1=os.path.join(self.basedir,self.source1,'std',self.var_name+'_'+str(self.level)+'_'+self.wildcard+'.nc')
            self.fileout1=os.path.join(self.basedir,self.source2,'std',self.var_name+'_'+str(self.level)+'_'+self.wildcard+'.nc')
        elif self.subdir=='processed':
            # one-off filein1 and fileout1 should have been passed as arguments
            pass
        else:
            raise UserWarning('Invalid subdir.')
        self.data_in=iris.load(self.filein1,self.name)
        if self.subdir=='std':
            # Get first time in input data
            x1=self.data_in[0].coord('time')[0]
            self.time1=x1.cell(0)[0]
            # Get last time in input data
            x1=self.data_in[-1].coord('time')[-1]
            self.time2=x1.cell(0)[0]
            self.time_units=x1.units
        if self.verbose:
            print(self)        

    def __repr__(self):
        return 'ModifySource({0.descriptor!r},verbose={0.verbose!r})'.format(self)

    def __str__(self):
        if self.verbose==2:
            ss=h1a+'ModifySource instance \n'+\
                'data_in: {0.data_in!s} \n'+\
                'source1: {0.source1!s} \n'+\
                'source2: {0.source2!s} \n'+\
                'frequency: {0.frequency!s} \n'+\
                'filein1: {0.filein1!s} \n'+h1b
            return(ss.format(self))
        else:
            return 'ModifySource instance'

    def f_time_average(self,method=3):
        """Time average data.

        Called from time_average.py.
        
        Used to change (reduce) time resolution of data, for subsequent
        data analysis.  Because the source attribute of a data set
        contains information on the time resolution, this effectively
        creates data with a different source.
        
        N.B. To calculate time mean statistics over a particular time
        domain, use the TimeDomStats class instead.
        
        N.B.  To "increase" time resolution of data, e.g., from weekly to
        daily, use the f_interpolate_time method in this class.
        """
        # Extract input data for current block of time
        time1,time2=block_times(self,verbose=self.verbose)
        if fnmatch.fnmatch(self.source1,'imerg???_sfc_30m') and fnmatch.fnmatch(self.source2,'imerg???_sfc_3h'):
            # Special case. Time averaging 30 minute IMERG data onto 3 hour TRMM time axis.
            # See comments below in method 3 for explanation
            tdc=datetime.timedelta(minutes=60)
            time1=time1-tdc
            time2=time2-tdc
            print('Modified time1,time2 for special case: {0!s}, {1!s}'.format(time1,time2))
        time_constraint=set_time_constraint(time1,time2,calendar=self.calendar,verbose=self.verbose)
        x1=self.data_in.extract(time_constraint)
        self.cube_in=x1.concatenate_cube()
        time_units=self.cube_in.coord('time').units
        
        if self.frequency in ['d','3h']:
            # Creating daily or 3h average data
            if method==1:
                print('f_time_average: Method 1')
                raise UserWarning('Do not use this method. Use method 3.')
                # Method 1 loops over days, creates a time constraint for
                # each day, extracts data for that day, then averages for that
                # day, then appends to a cube list, then finally merges the
                # cube list to create a cube of daily averaged data.  It is
                # very slow.
                #
                # Redundant code left here. Could remove.
                timedelta_day=datetime.timedelta(days=1)
                timedelta_minute=datetime.timedelta(seconds=60)
                timec1=time1
                # Create empty CubeList
                x10=iris.cube.CubeList([])
                while timec1<time2:
                    # Extract data over current day
                    timec2=timec1+timedelta_day-timedelta_minute
                    print(timec1,timec2)
                    time_constraintc=set_time_constraint(timec1,timec2,calendar=self.calendar,verbose=self.verbose)
                    x1=self.data_in.extract(time_constraintc)
                    x2=x1.concatenate_cube()
                    # Calculate daily mean
                    x3=x2.collapsed('time',iris.analysis.MEAN)
                    # Reset auxiliary time coordinate for current day at 00 UTC
                    timec_val=time_units.date2num(timec1)
                    timec_coord=iris.coords.DimCoord(timec_val,standard_name='time',units=time_units)
                    x3.remove_coord('time')
                    x3.add_aux_coord(timec_coord)
                    # Append current daily mean to cube list
                    x10.append(x3)
                    # Increment time
                    timec1+=timedelta_day
                x11=x10.merge_cube()
            elif method==2:
                print('f_time_average: Method 2')
                raise UserWarning('Do not use this method. Use method 3.')
                # Method 2 slices the input cube numpy array to get a
                # different numpy array for each time of day, then
                # adds them together, then divides by number of times
                # of day to get daily mean. Order(1000) faster than method 1
                # Method 2 depends on the time values being equally spaced,
                # and there being no missing data (in time)
                #
                # Redundant code left here. Could remove.
                #
                # Find time resolution of input data
                npd=find_npd(self.source1)
                # Check input data is well formed
                dim_coord_names=[xx.var_name for xx in self.cube_in.dim_coords]
                time_index=dim_coord_names.index('time')
                print('time_index: {0!s}'.format(time_index))
                if time_index!=0:
                    raise ToDoError('Code below only works if time is first dimension.')
                ntime=self.cube_in.shape[time_index]
                nday=int(ntime/npd)
                print('ntime,nday: {0!s}, {1!s}'.format(ntime,nday))
                if nday!=ntime/npd:
                    raise UserWarning('Input data is not integer number of days.')
                # Slice numpy array, one slice per time of day.
                # Then calculate daily means
                x1=self.cube_in.data
                x2=[x1[ii:ntime:npd,...] for ii in range(npd)]
                x3=x2[0]
                for xx in x2[1:]:
                    x3+=xx
                x4=x3/npd
                # Create new time axis for daily mean data (00 UTC each day)
                timec=time1
                time_vals=[]
                timedelta_day=datetime.timedelta(days=1)
                while timec<time2:
                    time_vals.append(time_units.date2num(timec))
                    timec+=timedelta_day
                time_coord=iris.coords.DimCoord(time_vals,standard_name='time',units=time_units)
                # Create new cube of daily mean
                x11=create_cube(x4,self.cube_in,new_axis=time_coord)
                cm=iris.coords.CellMethod('point','time',comments='daily mean from f_time_mean method 2')
                x11.add_cell_method(cm)
                print('Changing frequency attribute from {0!s} to {1!s}'.format(x11.attributes['frequency'],self.frequency))
                x11.attributes['frequency']=self.frequency
            elif method==3:
                print('f_time_average: Method 3')
                # Method 3 is identical to method 2 but has been generalised to
                # allow averaging to time resolution other than daily
                # This led to renaming of npd to nave, and nday to ntime2
                # and recoding how the new nave is found
                # and recoding creation of new time axis
                # Left method 2 code intact in case of error here.
                #
                # NB Time stamp convention for source2 is time stamp at
                # BEGINNING of interval that has been averaged over.
                # E.g., source1 is 'd' and source2 is '3h'
                # Will average over all 3-hourly data in a given calendar day
                # ie 00, 03, 06, 09, 12, 15, 18, 21 UTC
                # and write this with a time stamp of 00 UTC on that same
                # calendar day
                #
                # ---------------------------------------------------------------
                # However, there is one exception to the time stamp convention.
                # It is used to time average GPM IMERG 30-minute data onto the
                # TRMM 3-hour time axis, to effectively extend the TRMM data set
                # beyond Dec 2019, when the TRMM data ended.
                #
                # For this special case source1 must be 'imergXXX_sfc_30m' and 
                # source2 must be 'imergXXX_sfc_3h'
                #
                # NB Time stamp convention for source2 is time stamp at
                # approximate CENTRE of interval that has been averaged over.
                #
                # 30 minute IMERG data (source1) has times 0000, 0030, 0100, 0130,
                #  0200, 0230, ..., 2200, 2230, 2330 UTC
                # The target 3 hour time averaged data (source2) has times 0000,
                #  0300, 0600, 0900, 1200, 1500, 1800, 2100 UTC
                # E.g., to calculate a 0900 UTC value for source2 do an average
                #  of the 30 minute source1 values at  0800, 0830, 0900, 0930, 1000, 1030
                # Note that the actual central time of these input values is 0915 not 0930
                #  but this error is considered acceptable for typical uses of this data set
                # This means that special care must be taken for averaged values at the 
                #  beginning of each day 0000, and end of each day 2100
                # Average for 0000 is average of
                #  2300, 2330 of previous day, 0000, 0030, 0100, 0130.
                #  So at beginning of block need to read in from 2300, 2330 the previous day
                # Average for 2100 is average of
                #  2000, 2030, 2100, 2130, 2200, 2230
                #  So at end of block do not need 2300, 2330 in final day.
                # This explains offsetting of time1 and time2 in reading of block at beginning
                # of this function.
                # ------------------------------------------------------------
                #
                # Find time resolution of input and output data
                npd1=find_npd(self.source1)
                npd2=find_npd(self.source2)
                nave=int(npd1/npd2)
                print('npd1,npd2,nave: {0!s}, {1!s}, {2!s}'.format(npd1,npd2,nave))
                if nave!=npd1/npd2:
                    raise UserWarning('nave should be an integer.')
                # Check input data is well formed
                dim_coord_names=[xx.var_name for xx in self.cube_in.dim_coords]
                time_index=dim_coord_names.index('time')
                print('time_index: {0!s}'.format(time_index))
                if time_index!=0:
                    raise UserWarning('Time must be first dimension.')
                ntime=self.cube_in.shape[time_index]
                ntime2=int(ntime/nave)
                print('ntime,ntime2: {0!s}, {1!s}'.format(ntime,ntime2))
                if ntime2!=ntime/nave:
                    raise UserWarning('Input data is not integer number of source2 time interval.')
                # Slice numpy array, one slice per time of source2 time interval (day if 'd').
                # Then calculate source2 time interval means
                x1=self.cube_in.data
                x2=[x1[ii:ntime:nave,...] for ii in range(nave)]
                x3=x2[0]
                for xx in x2[1:]:
                    x3+=xx
                x4=x3/nave
                # Create new time axis for data on source2 time interval
                timec=time1
                if fnmatch.fnmatch(self.source1,'imerg???_sfc_30m') and fnmatch.fnmatch(self.source2,'imerg???_sfc_3h'):
                    # Special case. Averaging IMERG 30 min data to TRMM 3 hour time axis
                    # Reset starting time so new time axis is at approx CENTRE of input times and corresponds
                    # to TRMM 3 hour time axis
                    timec=time1+tdc
                time_vals=[]
                if self.frequency=='d':
                    timedelta_source2=datetime.timedelta(days=1)
                elif self.frequency=='3h':
                    timedelta_source2=datetime.timedelta(hours=3)
                else:
                    raise ToDoError('Code for other frequency.')
                while timec<time2:
                    time_vals.append(time_units.date2num(timec))
                    timec+=timedelta_source2
                time_coord=iris.coords.DimCoord(time_vals,standard_name='time',units=time_units)
                # Create new cube of data averaged onto source2 time interval
                x11=create_cube(x4,self.cube_in,new_axis=time_coord)
                cm=iris.coords.CellMethod('point','time',comments='daily mean from f_time_mean method 3')
                x11.add_cell_method(cm)
                print('Changing frequency attribute from {0!s} to {1!s}'.format(x11.attributes['frequency'],self.frequency))
                x11.attributes['frequency']=self.frequency
            else:
                raise UserWarning('Invalid method option.')
        else:
            raise ToDoError('The code above should work with frequencies other than d or 3h, but check.')
        # Convert units for selected data sources
        if self.source1 in ['trmm3b42v7_sfc_3h',] and self.source2 in ['trmm3b42v7_sfc_d',]:
            print("Converting TRMM precipitation from 3-hourly in 'mm hr-1' to daily mean in 'mm day-1'")
            x11.convert_units('mm day-1')
        self.cube_out=x11
        # Save time averaged cube
        if self.outfile_frequency=='year':
            x2=str(self.year)
        elif self.outfile_frequency=='month':
            x2=str(self.year)+str(self.month).zfill(2)
        else:
            raise UserWarning('outfile_frequency not recognised')
        self.fileout1=os.path.join(self.basedir,self.source2,'std',self.var_name+'_'+str(self.level)+'_'+x2+'.nc')
        iris.save(self.cube_out,self.fileout1)
        if self.archive:
            archive_file(self,self.fileout1)
        
    def f_interpolate_time(self):
        """Interpolate over time.

        Called from interpolate.py.
        
        N.B.  To decrease time resolution of data, e.g., from 3-hourly to
        daily, use the f_time_average method in this class.
        """
        # Ensure time interval to interpolate onto is within input time interval
        if self.time1_out<self.time1:
            self.time1_out=self.time1
            print('Correcting time1_out to be within bounds: {0.time1_out!s}'.format(self))
        if self.time2_out>self.time2:
            self.time2_out=self.time2
            print('Correcting time2_out to be within bounds: {0.time2_out!s}'.format(self))
        # Interpolate to daily data
        if self.frequency=='d':
            # Set times to interpolate onto
            self.time1_out_val=self.time_units.date2num(self.time1_out)
            self.time2_out_val=self.time_units.date2num(self.time2_out)
            if 'day' in self.time_units.name:
                time_diff=1
                self.sample_points=[('time',np.arange(self.time1_out_val,self.time2_out_val,time_diff))]
            else:
                raise ToDoError('Need code for time units that are not in days.')
            # Create cube (from cube list) of input data for interpolation
            # Times interval for this cube must completely contain the
            #   output interval (time_val1 to time_val2) for interpolation,
            #   if possible, otherwise there will be extrapolation at ends
            # Create a timedelta of default 25 days to account for this
            self.timedelta=datetime.timedelta(days=25)
            self.time1_in=self.time1_out-self.timedelta
            self.time2_in=self.time2_out+self.timedelta
            ss=h2a+'f_time_interpolate. \n'+\
                'timedelta: {0.timedelta!s} \n'+\
                'time1_in: {0.time1_in!s} \n'+\
                'time1_out: {0.time1_out!s} \n'+\
                'time2_out: {0.time2_out!s} \n'+\
                'time2_in: {0.time2_in!s} \n'+\
                'time1_out_val: {0.time1_out_val!s} \n'+\
                'time2_out_val: {0.time2_out_val!s} \n'+h2b
            if self.verbose==2:
                print(ss.format(self))
            time_constraint=set_time_constraint(time1_in,time2_in,calendar=self.calendar,verbose=self.verbose)
            x1=self.data_in.extract(time_constraint)
            self.cube_in=x1.concatenate_cube()
            # Interpolate in time
            self.cube_out=iris.cube.Cube.interpolate(self.cube_in,self.sample_points,scheme=iris.analysis.Linear())
            # Add cell method to describe linear interpolation
            cm=iris.coords.CellMethod('point','time',comments='linearly interpolated from weekly to daily time dimension')
            self.cube_out.add_cell_method(cm)
            # Save interpolated data
            fileout=replace_wildcard_with_time(self,self.fileout1)
            print('fileout: {0!s}'.format(fileout))
            iris.save(self.cube_out,fileout)
            if self.archive:
                archive_file(self,fileout)
        else:
            raise ToDoError('Need code for interpolation to other than daily data.')

    def f_zero_missing_values(self):
        """Set missing values to zero and save to new source.

        Called from zero_missing_values.py.
        """

        # Extract input data for current block of time
        time1,time2=block_times(self,verbose=self.verbose)
        time_constraint=set_time_constraint(time1,time2,calendar=self.calendar,verbose=self.verbose)
        x1=self.data_in.extract(time_constraint)
        self.cube_in=x1.concatenate_cube()
        # Set missing values to zero
        counter,self.cube_out=create_counter_from_mask(self.cube_in)
        counter_sum=counter.sum()
        npts=np.prod(counter.shape)
        fraction=counter_sum/npts
        print('counter_sum, npts, fraction: {0!s}, {1!s} {2!s}'.format(counter_sum,npts,fraction))
        # Add cell method to describe linear interpolation
        cm=iris.coords.CellMethod('point','time',comments='Masked values set to zero.')
        self.cube_out.add_cell_method(cm)
        # Save interpolated data
        fileout=replace_wildcard_with_time(self,self.fileout1)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.cube_out,fileout)
        if self.archive:
            archive_file(self,fileout)

    def f_get_target_grid(self):
        """Get grid from sample input data to use in later regridding.

        Option to restrict range of latitudes in grid, using self.latmin and self.latmax.

        Create attribute self.cube_grid.

        Called from regrid.py.
        """

        # Read sample data
        x1=iris.load(self.file_grid)
        self.cube_grid=x1.concatenate_cube()
        if not(self.latmin==self.latmax==False):
            print('# Restrict latitude range of grid.')
            atol=0.01
            print('latmin,latmax,atol: {0!s}, {1!s}, {2!s}'.format(self.latmin,self.latmax,atol))
            lat_constraint=iris.Constraint(latitude=lambda cell: self.latmin-atol <=cell<= self.latmax+atol)
            self.cube_grid=self.cube_grid.extract(lat_constraint)
            latcoord=self.cube_grid.coord('latitude')
            print('Latitude coord of new grid: {0!s}'.format(latcoord))
        self.cube_grid.coord(axis='x').guess_bounds()
        self.cube_grid.coord(axis='y').guess_bounds()
        if not(self.cube_grid.coord(axis='x').coord_system):
            # If data does not have a CoordSystem set, then give it one
            # as this is needed for regridding.
            # Assumed here that any data without a CoordSystem set will be on 
            # a latitude-longitude grid.
            self.new_cs=iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
            self.cube_grid.coord(axis='x').coord_system=self.new_cs
            self.cube_grid.coord(axis='y').coord_system=self.new_cs
        else:
            self.new_cs=self.cube_grid.coord(axis='x').coord_system

    def f_regrid(self):
        """Regrid data and save output.

        Create attribute self.cube_out.

        Called from regrid.py.
        """

        if self.subdir=='std':
            # Extract input data for current block of time
            time1,time2=block_times(self,verbose=self.verbose)
            time_constraint=set_time_constraint(time1,time2,calendar=self.calendar,verbose=self.verbose)
            x1=self.data_in.extract(time_constraint)
        elif self.subdir=='processed':
            x1=self.data_in
        else:
            raise UserWarning('Invalid subdir.')
        self.cube_in=x1.concatenate_cube()
        self.cube_in.coord(axis='x').guess_bounds()
        self.cube_in.coord(axis='y').guess_bounds()
        self.cube_in.coord(axis='x').coord_system=self.new_cs
        self.cube_in.coord(axis='y').coord_system=self.new_cs

        # Regrid
        self.cube_out=self.cube_in.regrid(self.cube_grid,iris.analysis.AreaWeighted())

        # Save regridded cube
        if self.subdir=='std':
            if self.outfile_frequency=='year':
                x2=str(self.year)
            elif self.outfile_frequency=='month':
                x2=str(self.year)+str(self.month).zfill(2)
            else:
                raise UserWarning('outfile_frequency not recognised')
            fileoutc=self.fileout1.replace(self.wildcard,x2)
        elif self.subdir=='processed':
            fileoutc=self.fileout1
        iris.save(self.cube_out,fileoutc)
        if self.archive:
            archive_file(self,fileoutc)

#==========================================================================

class SpatialSubset(object):
    
    """Create a SpatialSubset object (Hovmoller or 1-d time series).

    Called from spatial_subset.py

    Input is a N-dimensional cube, with one of these dimensions
    possibly being time, and the remaining dimensions being spatial
    coordinates.  Extract a spatial subset over one or more spatial
    dimensions and average over these dimensions.

    If all spatial dimensions are chosen (and subsetted and averaged
    over), and input had a time dimension, output is a 1-dimensional
    time series (area-averaged over a spatial box).

    If not all spatial dimensions are chosen, and input contains a
    time dimension, output contains a time dimension plus at least one
    spatial dimension, and is essentially a Hovmoller object.

    Selected attributes:

    self.data_in : iris cube list of all input data.

    self.data_subset_current : iris cube of subset of current data
    (e.g., year).

    self.data_subset : iris cube list of all subsetted data.

    self.file_data_in : path name for file(s) of input data.  If input
    data has a time dimension, contains a wild card * character, which
    will be replaced by, e.g., year numbers (if self.outfile_frequency
    is 'year').

    self.file_data_subset : path name for file(s) of subsetted data.
    If input data has a time dimension, contains a wild card *
    character, which will be replaced by, e.g., year numbers (if
    self.outfile_frequency is 'year').

    self.[band1_name,band1_val1,band1_val2] : Input cube is subset
    over dimension band1_name, and averaged between band1_val1 and
    band1_val2

    self.band2 : logical switch to enable subsetting and averaging
    over a second spatial dimension.

    self.[band2_name,band2_val1,band2_val2] : as for band1

    self.band3 : does not exist yet, but code could easily be extended
    to allow for subsetting over 3 spatial dimensions.
    
    History. 5 Jan 2018. This is a rename and generalisation of
    discontinued class Hovmoller.
    """

    def __init__(self,**descriptor):
        """Initialise from descriptor dictionary.

        Compulsory keywords: 'verbose','source','var_name','level'
        'basedir','subdir','filepre','band1_name','band1_val1','band1_val2','band2'.

        Optional keywords: 'band2_name','band2_val1','band2_val2'
        """
        self.__dict__.update(descriptor)
        self.descriptor=descriptor
        self.name=var_name2long_name[self.var_name]
        source_info(self)
        if self.subdir=='std':
            # Input data assumed to contain time dimension
            self.file_data_in=os.path.join(self.basedir,self.source,self.subdir,self.var_name+'_'+str(self.level)+self.filepre+'_'+self.wildcard+'.nc')
        elif self.subdir=='processed':
            # Input data does not contain a time dimension
            self.file_data_in=os.path.join(self.basedir,self.source,self.subdir,self.var_name+'_'+str(self.level)+self.filepre+'.nc')
        else:
            raise UserWarning('subdir is invalid.')
        self.strsubset='_ss_'+self.band1_name[:3]+'_'+str(self.band1_val1)+'_'+str(self.band1_val2)
        if self.band2:
            self.strsubset+='_'+self.band2_name[:3]+'_'+str(self.band2_val1)+'_'+str(self.band2_val2)
        if self.subdir=='std':
            self.file_data_subset=os.path.join(self.basedir,self.source,'processed',self.var_name+'_'+str(self.level)+self.filepre+self.strsubset+'_'+self.wildcard+'.nc')
        elif self.subdir=='processed':
            self.file_data_subset=os.path.join(self.basedir,self.source,'processed',self.var_name+'_'+str(self.level)+self.filepre+self.strsubset+'.nc')
        self.data_in=iris.load(self.file_data_in,self.name)
        if self.verbose:
            print(self)        

    def __repr__(self):
        return 'SpatialSubset({0.descriptor!r},verbose={0.verbose!r})'.format(self)

    def __str__(self):
        if self.verbose==2:
            ss=h1a+'SpatialSubset instance \n'+\
                'file_data_in: {0.file_data_in!s} \n'+\
                'data_in: {0.data_in!s} \n'+\
                'source: {0.source!s} \n'+\
                'band1_name: {0.band1_name!s} \n'+\
                'band1_val1: {0.band1_val1!s} \n'+\
                'band1_val2: {0.band1_val2!s} \n'
            if self.band2:
                ss+='band2_name: {0.band2_name!s} \n'+\
                     'band2_val1: {0.band2_val1!s} \n'+\
                     'band2_val2: {0.band2_val2!s} \n'
            ss+='file_data_subset: {0.file_data_subset!s} \n'+h1b
            return(ss.format(self))
        else:
            return 'SpatialSubset instance'

    def _constraint_(self,band_name,band_val1,band_val2,atol=1e-3):
        """Create constraint for band_name between band_val1 and band_val2. """
        if band_name=='latitude':
            band_constraint=iris.Constraint(latitude=lambda cell: band_val1-atol <=cell<= band_val2+atol)
        elif band_name=='longitude':
            band_constraint=iris.Constraint(longitude=lambda cell: band_val1-atol <=cell<= band_val2+atol)
        else:
            raise UserWarning('Invalid band_name.')
        return band_constraint

    def f_spatial_subset(self):
        """Create cube of spatially subsetted data and save.

        Create data_subset_current attribute.
        
        """
        # Extract input data for current time block (if appropriate), and 
        # for dimension band1_name, between band1_val1 and band1_val2
        # and optionally
        # for dimension band2_name, between band2_val1 and band2_val2
        if self.subdir=='std':
            self.time1,self.time2=block_times(self,verbose=self.verbose)
            time_constraint=set_time_constraint(self.time1,self.time2,calendar=self.calendar,verbose=self.verbose)
            xx1=self.data_in.extract(time_constraint)
        elif self.subdir=='processed':
            xx1=self.data_in
        band1_constraint=self._constraint_(self.band1_name,float(self.band1_val1),float(self.band1_val2))
        if self.band2:
            band2_constraint=self._constraint_(self.band2_name,float(self.band2_val1),float(self.band2_val2))            
            xx1=xx1.extract(band1_constraint & band2_constraint)
        else:
            xx1=xx1.extract(band1_constraint)
        ncubes=len(xx1)
        if ncubes!=1:
            raise UserWarning('Not a single cube. ncubes='+str(ncubes))
        xx1=xx1.concatenate_cube()
        self.cube=xx1
        print('self.cube.shape: {0!s}'.format(self.cube.shape))
        # Average over required spatial dimensions if necessary
        if self.band1_val1!=self.band1_val2:
            xx1=xx1.collapsed(self.band1_name,iris.analysis.MEAN)
            print('xx1.shape after averaging over band1: {0!s}'.format(xx1.shape))
            str1='Mean over {0.band1_name!s} {0.band1_val1!s} to {0.band1_val2!s}'.format(self)
            if self.band2 and self.band2_val1!=self.band2_val2:
                xx1=xx1.collapsed(self.band2_name,iris.analysis.MEAN)
                print('xx1.shape after averaging over band2: {0!s}'.format(xx1.shape))
                str1+=', and mean over {0.band2_name!s} {0.band2_val1!s} to {0.band2_val2!s}'.format(self)
            # Add a cell method to further describe averaging
            cm=iris.coords.CellMethod('mean',str1)
            if self.verbose:
                print(cm)
            xx1.add_cell_method(cm)
        self.data_subset_current=xx1
        # Save subsetted data
        if self.subdir=='std':
            fileout=replace_wildcard_with_time(self,self.file_data_subset)
        elif self.subdir=='processed':
            fileout=self.file_data_subset
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.data_subset_current,fileout)
        if self.archive:
            archive_file(self,fileout)

    def f_read_spatial_subset(self):
        """Read previously calculated spatially subsetted data.

        Create data_subset attribute.
        """
        self.data_subset=iris.load(self.file_data_subset,self.name)

#==========================================================================

class Wind(object):
    
    """Create a Wind object.

    Called from wind.py.

    Wrapper for windspharm VectorWind object to handle file i/o.

    Attributes (following VectorWind):

    self.uwnd : iris cube of u wind
    self.vwnd : iris cube of v wind
    self.psi : iris cube of streamfunction
    self.chi : iris cube of velocity potential
    self.vrt : iris cube of relative vorticity
    self.div : iris cube of divergence
    self.wndspd : iris cube of wind speed
    self.uwndchi : iris cube of eastward component of irrotational wind
    self.vwndchi : iris cube of northward component of irrotational wind
    self.duwnddx : iris cube of zonal derivative of eastward wind
    self.dvwnddy : iris cube of meridional derivative of northward wind

    self.flag_psi : Boolean flag to compute streamfunction
    or not.

    Similarly, self.flag_chi, self.flag_vrt, self.flag_div, etc.

    self.subdir : 'std', data has time dimension. Process for each time
                  'processed', data is time-invariant, e.g., time mean

    self.file_data : path name for file(s) of input data. if If input
    data has a time dimension, contains a wild card * character, which
    will be replaced by, e.g., year numbers (if self.outfile_frequency
    is 'year').  Contains a dummy string'VARNAME' to be replaced by
    'uwnd', 'vwnd', 'psi', etc.

    """

    def __init__(self,**descriptor):
        """Initialise from descriptor dictionary.

        Compulsory keywords: 'verbose','source','level',
        'basedir','subdir','filepre',
        'flag_psi','flag_chi','flag_vrt','flag_div',
        'flag_windspd', 'flag_uchi', 'flag_vchi'
        """
        self.__dict__.update(descriptor)
        self.descriptor=descriptor
        source_info(self)
        if self.subdir=='std':
            self.file_data=os.path.join(self.basedir,self.source,self.subdir,'VAR_NAME_'+str(self.level)+self.filepre+'_'+self.wildcard+'.nc')
        elif self.subdir=='processed':
            self.file_data=os.path.join(self.basedir,self.source,self.subdir,'VAR_NAME_'+str(self.level)+self.filepre+'.nc')
        else:
            raise UserWarning('subdir in invalid.')
        # uwnd
        self.var_name_uwnd='uwnd'
        self.name_uwnd=var_name2long_name[self.var_name_uwnd]
        self.file_data_uwnd=self.file_data.replace('VAR_NAME',self.var_name_uwnd)
        # vwnd
        self.var_name_vwnd='vwnd'
        self.name_vwnd=var_name2long_name[self.var_name_vwnd]
        self.file_data_vwnd=self.file_data.replace('VAR_NAME',self.var_name_vwnd)
        # psi
        if self.flag_psi:
            self.var_name_psi='psi'
            self.name_psi=var_name2long_name[self.var_name_psi]
            self.file_data_psi=self.file_data.replace('VAR_NAME',self.var_name_psi)
        # chi
        if self.flag_chi:
            self.var_name_chi='chi'
            self.name_chi=var_name2long_name[self.var_name_chi]
            self.file_data_chi=self.file_data.replace('VAR_NAME',self.var_name_chi)
        # vrt
        if self.flag_vrt:
            self.var_name_vrt='vrt'
            self.name_vrt=var_name2long_name[self.var_name_vrt]
            self.file_data_vrt=self.file_data.replace('VAR_NAME',self.var_name_vrt)
        # div
        if self.flag_div:
            self.var_name_div='div'
            self.name_div=var_name2long_name[self.var_name_div]
            self.file_data_div=self.file_data.replace('VAR_NAME',self.var_name_div)
        # wndspd
        if self.flag_wndspd:
            self.var_name_wndspd='wndspd'
            self.name_wndspd=var_name2long_name[self.var_name_wndspd]
            self.file_data_wndspd=self.file_data.replace('VAR_NAME',self.var_name_wndspd)
        # uwndchi
        if self.flag_uwndchi:
            self.var_name_uwndchi='uwndchi'
            self.name_uwndchi=var_name2long_name[self.var_name_uwndchi]
            self.file_data_uwndchi=self.file_data.replace('VAR_NAME',self.var_name_uwndchi)
        # vwndchi
        if self.flag_vwndchi:
            self.var_name_vwndchi='vwndchi'
            self.name_vwndchi=var_name2long_name[self.var_name_vwndchi]
            self.file_data_vwndchi=self.file_data.replace('VAR_NAME',self.var_name_vwndchi)
        # duwnddx
        if self.flag_duwnddx:
            self.var_name_duwnddx='duwnddx'
            self.name_duwnddx=var_name2long_name[self.var_name_duwnddx]
            self.file_data_duwnddx=self.file_data.replace('VAR_NAME',self.var_name_duwnddx)
        # dvwnddy
        if self.flag_dvwnddy:
            self.var_name_dvwnddy='dvwnddy'
            self.name_dvwnddy=var_name2long_name[self.var_name_dvwnddy]
            self.file_data_dvwnddy=self.file_data.replace('VAR_NAME',self.var_name_dvwnddy)
        #
        source_info(self)
        self.data_uwnd=iris.load(self.file_data_uwnd,self.name_uwnd)
        self.data_vwnd=iris.load(self.file_data_vwnd,self.name_vwnd)
        if self.verbose:
            print(self)        

    def __repr__(self):
        return 'Wind({0.descriptor!r},verbose={0.verbose!r})'.format(self)

    def __str__(self):
        if self.verbose==2:
            ss=h1a+'Wind instance \n'+\
                'source: {0.source!s} \n'+\
                'file_data: {0.file_data!s} \n'+\
                'file_data_uwnd: {0.file_data_uwnd!s} \n'+\
                'file_data_vwnd: {0.file_data_vwnd!s} \n'
            if self.flag_psi:
                ss+='file_data_psi: {0.file_data_psi!s} \n'
            if self.flag_chi:
                ss+='file_data_chi: {0.file_data_chi!s} \n'
            if self.flag_vrt:
                ss+='file_data_vrt: {0.file_data_vrt!s} \n'
            if self.flag_div:
                ss+='file_data_div: {0.file_data_div!s} \n'
            if self.flag_wndspd:
                ss+='file_data_wndspd: {0.file_data_wndspd!s} \n'
            if self.flag_uwndchi:
                ss+='file_data_uwndchi: {0.file_data_uwndchi!s} \n'
            if self.flag_vwndchi:
                ss+='file_data_vwndchi: {0.file_data_vwndchi!s} \n'
            if self.flag_duwnddx:
                ss+='file_data_duwnddx: {0.file_data_duwnddx!s} \n'
            if self.flag_dvwnddy:
                ss+='file_data_dvwnddy: {0.file_data_dvwnddy!s} \n'
            ss+=h1b
            return(ss.format(self))
        else:
            return 'Wind instance'

    def f_wind(self):
        """Calculate and save streamfunction etc.

        windspharm methods always return output with latitude running
        from north to south.  If input data runs south to north,
        switch the output data to also run from south to north, with
        lat_direction.
        """
        if self.subdir=='std':
            # Set current time range
            self.time1,self.time2=block_times(self,verbose=self.verbose)
            time_constraint=set_time_constraint(self.time1,self.time2,calendar=self.calendar,verbose=self.verbose)
            # Read uwnd and vwnd for current time range
            x1=self.data_uwnd.extract(time_constraint)
            x2=self.data_vwnd.extract(time_constraint)
        elif self.subdir=='processed':
            x1=self.data_uwnd
            x2=self.data_vwnd
        self.uwnd=x1.concatenate_cube()
        self.vwnd=x2.concatenate_cube()
        # Find value of south2north
        self.south2north=f_south2north(self.uwnd,verbose=True)
        # Create VectorWind instance
        self.ww=VectorWind(self.uwnd,self.vwnd)
        # Both psi and chi
        if self.flag_psi and self.flag_chi:
            self.psi,self.chi=self.ww.sfvp()
            self.psi.rename(self.name_psi)
            self.chi.rename(self.name_chi)
            self.psi.var_name=self.var_name_psi
            self.chi.var_name=self.var_name_chi
            if self.south2north:
                self.psi=lat_direction(self.psi,'s2n')
                self.chi=lat_direction(self.chi,'s2n')
            fileout1=self.file_data_psi
            fileout2=self.file_data_chi
            if self.subdir=='std':
                fileout1=replace_wildcard_with_time(self,fileout1)
                fileout2=replace_wildcard_with_time(self,fileout2)
            print('fileout1: {0!s}'.format(fileout1))
            print('fileout2: {0!s}'.format(fileout2))
            iris.save(self.psi,fileout1)
            iris.save(self.chi,fileout2)
            if self.archive:
                archive_file(self,fileout1)
                archive_file(self,fileout2)
        # or psi only
        elif self.flag_psi:
            self.psi=self.ww.streamfunction()
            self.psi.rename(self.name_psi)
            self.psi.var_name=self.var_name_psi
            if self.south2north:
                self.psi=lat_direction(self.psi,'s2n')
            fileout=self.file_data_psi
            if self.subdir=='std':
                fileout=replace_wildcard_with_time(self,fileout)
            print('fileout: {0!s}'.format(fileout))
            iris.save(self.psi,fileout)
            if self.archive:
                archive_file(self,fileout)
        # or chi only
        elif self.flag_chi:
            self.chi=self.ww.velocitypotential()
            self.chi.rename(self.name_chi)
            self.chi.var_name=self.var_name_chi
            if self.south2north:
                self.chi=lat_direction(self.chi,'s2n')
            fileout=self.file_data_chi
            if self.subdir=='std':
                fileout=replace_wildcard_with_time(self,fileout)
            print('fileout: {0!s}'.format(fileout))
            iris.save(self.chi,fileout)
            if self.archive:
                archive_file(self,fileout)
        # Both vrt and div
        if self.flag_vrt and self.flag_div:
            self.vrt,self.div=self.ww.vrtdiv()
            self.vrt.rename(self.name_vrt)
            self.vrt.rename(self.name_vrt)
            self.vrt.var_name=self.var_name_vrt
            self.div.var_name=self.var_name_div
            if self.south2north:
                self.vrt=lat_direction(self.vrt,'s2n')
                self.div=lat_direction(self.div,'s2n')
            fileout1=self.file_data_vrt
            fileout2=self.file_data_div
            if self.subdir=='std':
                fileout1=replace_wildcard_with_time(self,fileout1)
                fileout2=replace_wildcard_with_time(self,fileout2)
            print('fileout1: {0!s}'.format(fileout1))
            print('fileout2: {0!s}'.format(fileout2))
            iris.save(self.vrt,fileout1)
            iris.save(self.div,fileout2)
            if self.archive:
                archive_file(self,fileout1)
                archive_file(self,fileout2)
        # or vrt only
        elif self.flag_vrt:
            self.vrt=self.ww.vorticity()
            self.vrt.rename(self.name_vrt)
            self.vrt.var_name=self.var_name_vrt
            if self.south2north:
                self.vrt=lat_direction(self.vrt,'s2n')
            fileout=self.file_data_vrt
            if self.subdir=='std':
                fileout=replace_wildcard_with_time(self,fileout)
            print('fileout: {0!s}'.format(fileout))
            iris.save(self.vrt,fileout)
            if self.archive:
                archive_file(self,fileout)
        # or div only
        elif self.flag_div:
            self.div=self.ww.divergence()
            self.div.rename(self.name_div)
            self.div.var_name=self.var_name_div
            if self.south2north:
                self.div=lat_direction(self.div,'s2n')
            fileout=self.file_data_div
            if self.subdir=='std':
                fileout=replace_wildcard_with_time(self,fileout)
            print('fileout: {0!s}'.format(fileout))
            iris.save(self.div,fileout)
            if self.archive:
                archive_file(self,fileout)
        # wndspd
        if self.flag_wndspd:
            self.wndspd=self.ww.magnitude()
            self.wndspd.rename(self.name_wndspd)
            self.wndspd.var_name=self.var_name_wndspd
            if self.south2north:
                self.wndspd=lat_direction(self.wndspd,'s2n')
            fileout=self.file_data_wndspd
            if self.subdir=='std':
                fileout=replace_wildcard_with_time(self,fileout)
            print('fileout: {0!s}'.format(fileout))
            iris.save(self.wndspd,fileout)
            if self.archive:
                archive_file(self,fileout)
        # uwndchi and vwndchi
        if self.flag_uwndchi and self.flag_vwndchi:
            self.uwndchi,self.vwndchi=self.ww.irrotationalcomponent()
            self.uwndchi.rename(self.name_uwndchi)
            self.vwndchi.rename(self.name_vwndchi)
            self.uwndchi.var_name=self.var_name_uwndchi
            self.vwndchi.var_name=self.var_name_vwndchi
            if self.south2north:
                self.uwndchi=lat_direction(self.uwndchi,'s2n')
                self.vwndchi=lat_direction(self.vwndchi,'s2n')
            fileout1=self.file_data_uwndchi
            fileout2=self.file_data_vwndchi
            if self.subdir=='std':
                fileout1=replace_wildcard_with_time(self,fileout1)
                fileout2=replace_wildcard_with_time(self,fileout2)
            print('fileout1: {0!s}'.format(fileout1))
            print('fileout2: {0!s}'.format(fileout2))
            iris.save(self.uwndchi,fileout1)
            iris.save(self.vwndchi,fileout2)
            if self.archive:
                archive_file(self,fileout1)
                archive_file(self,fileout2)
        # duwnddx
        if self.flag_duwnddx:
            # Create dummy cube of zeros for vwnd
            zeros=np.zeros(self.uwnd.data.shape)
            zeros=create_cube(zeros,self.uwnd,new_var_name='dummy')
            # Create dummy VectorWind instance
            ww1=VectorWind(self.uwnd,zeros)
            # Calculate duwnddx as usual
            self.duwnddx=ww1.divergence()
            self.duwnddx.rename(self.name_duwnddx)
            self.duwnddx.var_name=self.var_name_duwnddx
            if self.south2north:
                self.duwnddx=lat_direction(self.duwnddx,'s2n')
            fileout=self.file_data_duwnddx
            if self.subdir=='std':
                fileout=replace_wildcard_with_time(self,fileout)
            print('fileout: {0!s}'.format(fileout))
            iris.save(self.duwnddx,fileout)
            if self.archive:
                archive_file(self,fileout)
        # dvwnddy
        if self.flag_dvwnddy:
            # Create dummy cube of zeros for uwnd
            zeros=np.zeros(self.uwnd.data.shape)
            zeros=create_cube(zeros,self.uwnd,new_var_name='dummy')
            # Create dummy VectorWind instance
            ww1=VectorWind(zeros,self.vwnd)
            # Calculate dvwnddy as usual
            self.dvwnddy=ww1.divergence()
            self.dvwnddy.rename(self.name_dvwnddy)
            self.dvwnddy.var_name=self.var_name_dvwnddy
            if self.south2north:
                self.dvwnddy=lat_direction(self.dvwnddy,'s2n')
            fileout=self.file_data_dvwnddy
            if self.subdir=='std':
                fileout=replace_wildcard_with_time(self,fileout)
            print('fileout: {0!s}'.format(fileout))
            iris.save(self.dvwnddy,fileout)
            if self.archive:
                archive_file(self,fileout)

#==========================================================================

class CombineLevels(object):
    
    """Create a CombineLevels object.

    Called from combine_levels.py.

    Attributes:

    self.subdir : 'std', data has time dimension. Process for each time
                  'processed', data is time-invariant, e.g., time mean

    self.file_data : path name for file(s) of input data. If input
    data has a time dimension, contains a wild card * character, which
    will be replaced by, e.g., year numbers (if self.outfile_frequency
    is 'year').

    """

    def __init__(self,**descriptor):
        """Initialise from descriptor dictionary.

        Compulsory keywords: 'verbose','source','levels',
        'basedir','subdir','filepre','var_name'
        """
        self.__dict__.update(descriptor)
        self.descriptor=descriptor
        source_info(self)
        if self.subdir=='std':
            self.file_data=os.path.join(self.basedir,self.source,self.subdir,self.var_name+'LEVEL'+self.filepre+'_'+self.wildcard+'.nc')
        elif self.subdir=='processed':
            self.file_data=os.path.join(self.basedir,self.source,self.subdir,self.var_name+'_LEVEL'+self.filepre+'.nc')
        else:
            raise UserWarning('subdir in invalid.')
        self.name=var_name2long_name[self.var_name]
        #
        self.data_in={}
        for levelc in self.levels:
            x1=self.file_data.replace('LEVEL',str(levelc))
            print(x1)
            self.data_in[self.var_name+'_'+str(levelc)]=iris.load(x1,self.name)
        if self.verbose:
            print(self)        

    def __repr__(self):
        return 'CombineLevels({0.descriptor!r},verbose={0.verbose!r})'.format(self)

    def __str__(self):
        if self.verbose==2:
            ss=h1a+'CombineLevels instance \n'+\
                'source: {0.source!s} \n'+\
                'file_data: {0.file_data!s} \n'
            ss+=h1b
            return(ss.format(self))
        else:
            return 'CombineLevels instance'

    def f_combine_levels(self):
        """Combine multi-level data into single cube.

        """
        # Read data at all levels and combine levels into single cube
        x2=iris.cube.CubeList([])
        if self.subdir=='std':
            self.time1,self.time2=block_times(self,verbose=self.verbose)
            time_constraint=set_time_constraint(self.time1,self.time2,calendar=self.calendar,verbose=self.verbose)
        for levelc in self.levels:
            self.level=levelc
            if self.subdir=='std':
                x1=self.data_in[self.var_name+'_'+str(self.level)].extract(time_constraint)
            else:
                x1=self.data_in[self.var_name+'_'+str(self.level)]
            x1=x1.concatenate_cube()
            levelcoord=iris.coords.AuxCoord(self.level,var_name='level')
            x1.add_aux_coord(levelcoord)
            x2.append(x1)
        # 14 May 2020, noted in f_combine_latitudes that order of dimensions is not conventional
        raise ToDoError('Check that dimensions are in normal order tzyx. Also that data is float32 not float64. See f_combine_latitudes for code and copy here.')
        self.data_all=x2.merge_cube()
        # Save multi-level data cube
        fileout1=self.file_data.replace('LEVEL','all')
        print('fileout1: {0!s}'.format(fileout1))
        iris.save(self.data_all,fileout1)
        if self.archive:
            archive_file(self,fileout1)

#==========================================================================

class CombineLatitudes(object):
    
    """Create a CombineLatitudes object.

    Called from combine_latitudes.py.

    An ad hoc class. Use when filtering data for equatorial
    waves. First use wheelerkiladis.py to: take
    time-latitude-longitude data split and stored by time (e.g.,
    yearly or monthly files), then split this into individual
    latitudes and recombine the time axis to get a single
    time-longitude Hovmollers at each latitude; filter these for a
    particular equatorial wave (e.g., equatorial Kelvin waves).

    Now use this class to recombine the filtered data. Recombine the
    latitudes, and resplit for storage by time (into yearly or monthly
    files again). Hence, this process creates a new "source" of data.

    Attributes:

    self.source1 : source of original data, e.g.,
    'erainterim_plev_6h'. The output from wheelerkiladis.py is stored
    in the 'processed' directory of this source.

    self.source2: source for the recombined data, e.g.,
    'erainterimEK1_plev_6h'.

    self.file_hovWKfilt : path name for files of input data. Contains
    the string 'LATITUDE', which will be replaced by individual
    latitudes.

    self.file_data_out : path name for file(s) of output
    data. Contains a wild card * character, which will be replaced by,
    e.g., year numbers (if self.outfile_frequency is 'year').

    """

    def __init__(self,**descriptor):
        """Initialise from descriptor dictionary.

        Compulsory keywords: 'verbose','source1', 'source2', 'level',
        'latitudes', 'basedir','filepre','var_name', 
        'wave_type, ''time1', 'time2'
        """
        self.__dict__.update(descriptor)
        self.descriptor=descriptor
        self.source=self.source1
        source_info(self)

        #ss='_lat_'+str(self.lat1)+'_'+str(self.lat2)+'_'+str(self.time1)[:10]+'_'+str(self.time2)[:10]
        self.file_hovWKfilt=os.path.join(self.basedir,self.source1,'processed',self.var_name+'_'+str(self.level)+self.filepre+'_hovWKfilt'+self.wave_type+'_lat_LATITUDE_'+str(self.time1)[:10]+'_'+str(self.time2)[:10]+'.nc')
        self.file_data_out=os.path.join(self.basedir,self.source2,'std',self.var_name+'_'+str(self.level)+self.filepre+'_'+self.wildcard+'.nc')
        self.name=var_name2long_name[self.var_name]
        #
        self.data_in={}
        for latitudec in self.latitudes:
            x1=self.file_hovWKfilt.replace('LATITUDE',str(latitudec)+'_'+str(latitudec))
            print(x1)
            self.data_in[self.var_name+'_'+str(latitudec)]=iris.load(x1,self.name)
        if self.verbose:
            print(self)        

    def __repr__(self):
        return 'CombineLatitudes({0.descriptor!r},verbose={0.verbose!r})'.format(self)

    def __str__(self):
        if self.verbose==2:
            ss=h1a+'CombineLatitudes instance \n'+\
                'latitudes: {0.latitudes!s} \n'+\
                'source1: {0.source1!s} \n'+\
                'source2: {0.source2!s} \n'+\
                'file_hovWKfilt: {0.file_hovWKfilt!s} \n'+\
                'file_data_out: {0.file_data_out!s} \n'
            ss+=h1b
            return(ss.format(self))
        else:
            return 'CombineLatitudes instance'

    def f_combine_latitudes(self):
        """Combine multi-latitude data into single cube.

        """
        # Read data at all latitudes and combine latitudes into single cube
        x2=iris.cube.CubeList([])
        self.time1,self.time2=block_times(self,verbose=self.verbose)
        time_constraint=set_time_constraint(self.time1,self.time2,calendar=self.calendar,verbose=self.verbose)
        for latitudec in self.latitudes:
            print('latitudec: {0!s}'.format(latitudec))
            self.latitude=latitudec
            x1=self.data_in[self.var_name+'_'+str(self.latitude)].extract(time_constraint)
            x1=x1.concatenate_cube()
            latitudecoord=iris.coords.AuxCoord(self.latitude,var_name='latitude')
            x1.add_aux_coord(latitudecoord)
            x2.append(x1)
        x2=x2.merge_cube()
        print('x2 before: {0!s}'.format(x2))
        #
        for xx in x2.aux_coords:
            if xx.var_name=='latitude_0':
                # Remove extraneous latitude auxilliary coordinate
                x2.remove_coord('latitude_0')
        #
        dc=x2.dim_coords
        dclen=len(dc)
        dcname0=dc[0].name()
        dcname1=dc[1].name()
        dcname2=dc[2].name()
        print('dclen,dcname0,dcname1,dcname2: {0!s}, {1!s}, {2!s}, {3!s}'.format(dclen,dcname0,dcname1,dcname2))
        if (dclen==3 and dcname0=='time' and dcname1=='latitude' and dcname2=='longitude'):
            # Dimensions are in order [time,latitude,longitude]
            pass
        if (dclen==3 and dcname0=='latitude' and dcname1=='time' and dcname2=='longitude'):
            # Dimensions are in order [latitude,time,longitude]
            # Reorder to [time,latitude,longitude]
            x2.transpose([1,0,2])
        else:
            raise ToDoError('Code up for dimension order other than time,lat,lon.')
        #
        print('dtype: {0!s}'.format(x2.data.dtype))
        if x2.data.dtype==np.float64:
            # Convert to float32 (single precision) as input is unnecessarily float62
            x2.data=np.float32(x2.data)
            print('dtype: {0!s}'.format(x2.data.dtype))
        self.data_all=x2
        print('x2 after: {0!s}'.format(x2))
        # Save multi-latitude data cube
        if self.outfile_frequency=='year':
            fileout1=self.file_data_out.replace(self.wildcard,str(self.year))
        elif self.outfile_frequency=='month':
            fileout1=self.file_data_out.replace(self.wildcard,str(self.year)+str(self.month).zfill(2))
        else:
            raise UserWarning('Invalid outfile_frequency.')
        print('fileout1: {0!s}'.format(fileout1))
        iris.save(self.data_all,fileout1)
        if self.archive:
            archive_file(self,fileout1)

#==========================================================================

class AnnualCycle(object):
    
    """Calculate and subtract annual cycle.

    Called from anncycle.py

    Assumes input data has equally spaced time intervals.

    Attributes:

    self.frequency : string to denote frequency of input (and output)
    data, e.g., 'd' for daily data.

    self.year1 and self.year2 : integers.  The annual cycle will be
    calculated using input data from 0000:00 UTC 1 Jan self.year1 to
    2359:59 31 Dec self.year2.

    self.data_in : iris cube list of all input data.

    self.data_anncycle_raw : iris cube of 'raw' annual cycle.  The
    data for e.g. 5 Jan is a simple mean of all the input data for 5
    Jan.

    self.kmin : integer >=1. Typical value is 5. Minimum number of
    valid input data values needed at a particular grid point on a
    particular day of the year, for the raw annual cycle to be
    calculated at that grid point and that day of year. If less than
    kmin data values are available, the annual cycle at that grid
    point and that day of year is set to the mean value of the whole
    global annual cycle. If only a small fraction (see self.frac_crit)
    of the raw annual cycle values at a grid point are set to this
    mean value, this will have an insignificant effect on the smoothed
    annual cycle calculated in f_anncycle_smooth().

    self.detrend : Boolean. If true, calculate the temporal linear
    trend over the whole period and remove this along with the mean
    and annual cycle.

    self.frac_crit: float. Currently hard set to 0.01. If greater than
    self.frac_crit of grid points have missing data in the raw annual
    cycle, raise UserWarning.  If this happens, examine carefully and
    recode accordingly.

    self.data_anncycle_smooth : iris cube of smoothed annual cycle.
    Calculated from mean plus first self.nharm annual harmonics of
    self.data_anncycle_raw.

    self.nharm : integer number of annual harmonics to retain for
    smoothed annual cycle.

    self.data_anncycle_rm : iris cube list of anomaly data with annual
    cycle subtracted.

    self.anncycle_source : Usually False. Creation of raw and smoothed
    annual cycle is currently only available for daily data (frequency
    of 'd'). If using daily data throughout, set anncycle_source to
    False. However, higher frequency (e.g., '6h') data can be
    processed. The raw and smoothed annual cycles must first be
    calculated separately from daily data (recommended to set detrend
    to False), using f_anncycle_raw and f_anncycle_smooth (with
    anncycle_source set to False). Then, run anncycle.py again, with
    SOURCE set to the higher frequency source (e.g.,
    'erainterim_plev_6h'), and anncycle_source set to the daily
    resolution source, e.g.,
    'erainterim_plev_d'. f_read_anncycle_smooth will be directed to
    read in the previously calculated smoothed annual cycle of daily
    data. This smoothed annual cycle of daily data is then expanded to
    the higher frequency using f_expand_anncycle_smooth, and can then
    be subtracted from the higher frequency e.g., 6h, input data,
    using f_subtract_anncycle.
    
    self.file_data_in : path name for file(s) of input data.  Contains a
    wild card ???? character, which will be replaced by, e.g., year
    numbers (if self.outfile_frequency is 'year').
    
    self.file_anncycle_mm : path name for file of linear trend
    gradient mm.
    
    self.file_anncycle_cc : path name for file of linear trend
    intercept cc.
    
    self.file_anncycle_raw : path name for file of raw annual cycle.
    
    self.file_anncycle_smooth : path name for file of smoothed annual
    cycle.
    
    self.file_anncycle_rm : path name for file(s) of data with
    smoothed annual cycle subtracted.  Contains a wild card ????
    character, which will be replaced by, e.g., year numbers (if
    self.outfile_frequency is 'year').

    self.time1: If False, calculate anomalies from annual cycle from
    beginning of data set.  If datetime.datetime object, calculate
    anomalies from annual cycle from this time.
    
    self.time2: If False, calculate anomalies from annual cycle up to
    end of data set.  If datetime.datetime object, calculate anomalies
    from annual cycle up to this time.

    """

    def __init__(self,**descriptor):
        """Initialise from descriptor dictionary.

        Compulsory keywords: 'verbose','source','var_name','level',
        'basedir','year1','year2','time1','time2','nharm','detrend', 
        'anncycle_source'.
        """
        self.__dict__.update(descriptor)
        self.descriptor=descriptor
        self.name=var_name2long_name[self.var_name]
        source_info(self)
        if self.detrend and self.anncycle_source:
            raise UserWarning('Stop. Have not considered impact of detrending on expanding annual cycle to higher frequency.')
        self.frac_crit=0.01
        self.file_data_in=os.path.join(self.basedir,self.source,'std',
              self.var_name+'_'+str(self.level)+'_'+self.wildcard+'.nc')
        self.file_anncycle_mm=os.path.join(self.basedir,self.source,
              'processed',self.var_name+'_'+str(self.level)+'_ac_mm_'+\
              str(self.year1)+'_'+str(self.year2)+'.nc')
        self.file_anncycle_cc=os.path.join(self.basedir,self.source,
              'processed',self.var_name+'_'+str(self.level)+'_ac_cc_'+\
              str(self.year1)+'_'+str(self.year2)+'.nc')
        self.file_anncycle_raw=os.path.join(self.basedir,self.source,
              'processed',self.var_name+'_'+str(self.level)+'_ac_raw_'+\
              str(self.year1)+'_'+str(self.year2)+'.nc')
        xx1=self.source
        if self.anncycle_source:
            xx1=self.anncycle_source
        self.file_anncycle_smooth=os.path.join(self.basedir,xx1,
              'processed',self.var_name+'_'+str(self.level)+'_ac_smooth_'+\
              str(self.year1)+'_'+str(self.year2)+'.nc')
        if self.calendar=='gregorian':
            self.file_anncycle_smooth_leap=self.file_anncycle_smooth.replace('.','_leap.')
        #
        if self.anncycle_source:
            self.file_expanded_anncycle_smooth=os.path.join(self.basedir,self.source,
                  'processed',self.var_name+'_'+str(self.level)+'_ac_smooth_expanded_'+\
                  str(self.year1)+'_'+str(self.year2)+'.nc')
            if self.calendar=='gregorian':
                self.file_expanded_anncycle_smooth_leap=self.file_expanded_anncycle_smooth.replace('.','_leap.')
        #
        self.file_anncycle_rm=os.path.join(self.basedir,self.source,'std',
              self.var_name+'_'+str(self.level)+'_rac_'+self.wildcard+'.nc')
        self.data_in=iris.load(self.file_data_in,self.name)
        self.file_individual_annual_means=os.path.join(self.basedir,self.source,'processed',
              self.var_name+'_'+str(self.level)+'_ann????.nc')
        if not(isinstance(self.kmin,int) and self.kmin>=1):
            raise UserWarning('kmin must be an integer greater or equal to 1')
        # yearac is dummy year to be used for time coordinate of annual cycle
        # Initially, was set to 1 (i.e., AD year 1). However, after a module
        # update in Mar 2019, this threw up errors in date comparisons.
        # Tracked to Gregorian calendar dates before 1582.
        # Solution, set yearac to a  year later than 1582. 
        # Choose 2101 (do not want a leap year)
        # Similarly with yearacleap, the dummy year for time coordinate of
        #  annual cycle in leap years
        #self.yearac=1; self.yearacleap=4
        self.yearac=2101; self.yearacleap=self.yearac+3
        # Get first time in input data if self.time1 not externally set
        if not self.time1:
            x1=self.data_in[0].coord('time')[0]
            x2=x1.cell(0)[0]
            if self.calendar=='gregorian':
                self.time1=create_datetime(x2)
            elif self.calendar=='360_day':
                self.time1=create_Datetime360Day(x2)
            else:
                raise UserWarning('Invalid calendar.')
        # Get last time in input data if self.time2 not externally set
        if not self.time2:
            x1=self.data_in[-1].coord('time')[-1]
            x2=x1.cell(0)[0]
            if self.calendar=='gregorian':
                self.time2=create_datetime(x2)
            elif self.calendar=='360_day':
                self.time2=create_Datetime360Day(x2)
        if self.verbose:
            print(self)        

    def __repr__(self):
        return 'AnnualCycle({0.descriptor!r},verbose={0.verbose!r})'.format(self)

    def __str__(self):
        if self.verbose==2:
            ss=h1a+'AnnualCycle instance \n'+\
                'year1: {0.year1!s} \n'+\
                'year2: {0.year2!s} \n'+\
                'time1: {0.time1!s} \n'+\
                'time2: {0.time2!s} \n'+\
                'nharm: {0.nharm!s} \n'+\
                'kmin: {0.kmin!s} \n'+\
                'detrend: {0.detrend!s} \n'+\
                'frac_crit: {0.frac_crit!s} \n'+\
                'file_data_in: {0.file_data_in!s} \n'+\
                'data_in: {0.data_in!s} \n'+\
                'file_anncycle_raw: {0.file_anncycle_raw!s} \n'+\
                'file_anncycle_smooth: {0.file_anncycle_smooth!s} \n'
            if self.calendar=='gregorian':
                ss+='file_anncycle_smooth_leap: {0.file_anncycle_smooth_leap!s} \n'
            ss+='file_anncycle_rm: {0.file_anncycle_rm!s} \n'+h1b
            return(ss.format(self))
        else:
            return 'AnnualCycle instance'

    def f_trend(self):
        """Calculate and save temporal linear trend over whole time period.

        Create mm and cc attributes.

        N.B., calculcation of the linear trend at each grid point of
        (typically) decades of daily (or sub-daily) data would either
        require unfeasible memory (to load in the entire data set)
        and/or would be computationally inefficient (by e.g., loop
        over longitude and latitude, and reading in the full time
        series for each grid point).

        Hence, a pragmatic approach is taken here.

        First, individual annual means are calculated, for each year
        between self.year1 and self.year2. Note that this must be done
        independently of (and before) running anncycle.py, using
        mean.py. These should be calculated using time domains with
        id's of the form annYYYY, e.g., ann1979.

        Hence, we then have a manageable number (typically 10-50) of
        individual annual mean fields (maps). These are read into
        memory, combined into a single cube, and the linear trend
        calculated at each grid point.

        The temporal linear trend is stored as maps self.mm and
        self.cc, to be used in an equation of the form yy = mm * time
        + cc. The trend time axis here is the same as in the input
        data, i.e., time=0 corresponds to the base time in the time
        coordinate of the input data, The constant cc corresponds to
        the value of the trend line at time=0.

        Later, in f_subtract_anncycle, the value yy is calculated and
        subtracted first, then the annual harmonics are
        subtracted. Note that the trend will include the time mean, so
        if the trend is subtracted, the time mean is not subtracted
        again when the annual harmonics are subtracted.


        """
        print('Read pre-calculated individual annual means.')
        x1=iris.cube.CubeList([])
        for yearc in range(self.year1,self.year2+1):
            filec=self.file_individual_annual_means.replace('????',str(yearc).zfill(4))
            print('filec: {0!s}'.format(filec))
            x2=iris.load(filec)
            x3=x2.concatenate_cube()
            # Replace time coordinate (to get rid of bounds, that prevent concatenation later)
            # Remove cell methods for same reason
            time_coord=x3.coord('time')
            time_coord2=iris.coords.AuxCoord(time_coord.points[0],standard_name='time',var_name='time',units=time_coord.units)
            x3.remove_coord('time')
            x3.add_aux_coord(time_coord2)
            x3.cell_methods=()
            x1.append(x3)
        self.ann_means=concatenate_cube(x1)
        print('Calculate linear temporal trend.')
        xx=self.ann_means.coord('time').points
        yy=self.ann_means.data
        xxshape=xx.shape
        yyshape=yy.shape
        print('xxshape, yyshape: {0!s}, {1!s}'.format(xxshape,yyshape))
        if self.ann_means.coords()[0].name()!='time':
            raise UserWarning('First dimension must be time.')
        ntime=yyshape[0]
        self.ngrid=np.prod(yyshape[1:])
        print('ntime, ngrid: {0!s}, {1!s}'.format(ntime,self.ngrid))
        yy2shape=(ntime,self.ngrid)
        yy2=yy.data.reshape(yy2shape)
        print('yy2shape: {0!s}'.format(yy2shape))
        zz=np.polyfit(xx,yy2,1)
        mm=zz[0]
        cc=zz[1]
        mm=mm.reshape(yyshape[1:])
        cc=cc.reshape(yyshape[1:])
        print('mm.shape, cc.shape {0!s}, {1!s}'.format(mm.shape,cc.shape))
        mm=create_cube(mm,x3,new_var_name='ltt')
        mm.units=str(self.ann_means.units)+' ('+self.frequency+')-1'
        mm.remove_coord('time')
        comments=str(time_coord.units)
        cm=iris.coords.CellMethod('mean','time',comments=comments)
        mm.add_cell_method(cm)
        self.mm=mm
        cc=create_cube(cc,x3)
        cc.remove_coord('time')
        self.cc=cc
        # Save linear trend
        iris.save(self.mm,self.file_anncycle_mm)
        iris.save(self.cc,self.file_anncycle_cc)
        if self.archive:
            archive_file(self,self.file_anncycle_mm)
            archive_file(self,self.file_anncycle_cc)

    def f_read_trend(self):
        """Read previously created linear trend.

        Create mm and cc attributes.
        
        """
        self.mm=iris.load_cube(self.file_anncycle_mm,'linear_temporal_trend')
        self.cc=iris.load_cube(self.file_anncycle_cc,self.name)

    def f_anncycle_raw(self):
        """Create raw annual cycle and write to file.

        The raw annual cycle consists of a value for each time (e.g.,
        day) of the year by averaging over all years, e.g., creates a
        value for 1 Jan by averaging over all the 1 Jan days from all
        the year.
        
        This will be assigned to year 1, i.e. AD 1.  NB year 1 is not
        a leap year.

        Ignores 29 Feb in leap years. Hence, final annual cycle has
        exactly 365 days of data.

        Accounts for missing data. See documentation of self.kmin and
        self.frac_crit in class docstring for explanation.

        Create data_anncycle_raw attribute.
        
        Note for future coding: First attempt at this, code was slow
        as it looped over each day of the year, and reads in all data
        (e.g., 1 Jan's) for that day and averages over them.

        """
        if self.frequency!='d':
            raise ToDoError('Annual cycle presently only coded up for daily data.')
        if self.calendar=='gregorian':
            # First year
            yearc=self.year1
            # Extract 1 Jan to 28 Feb of current year
            time1=datetime.datetime(yearc,1,1)
            time2=datetime.datetime(yearc,2,28)
            print(time1,time2)
            time_constraint=iris.Constraint(time=lambda cell: time1<=cell<=time2)
            xx1=self.data_in.extract(time_constraint)
            janfeb=xx1.concatenate_cube()
            janfeb_counter,janfeb=create_counter_from_mask(janfeb,verbose=self.verbose)
            janfeb.remove_coord('time')
            # Extract 1 Mar to 31 Dec of current year
            time1=datetime.datetime(yearc,3,1)
            time2=datetime.datetime(yearc,12,31)
            print(time1,time2)
            time_constraint=iris.Constraint(time=lambda cell: time1<=cell<=time2)
            xx1=self.data_in.extract(time_constraint)
            mardec=xx1.concatenate_cube()
            mardec_counter,mardec=create_counter_from_mask(mardec,verbose=self.verbose)
            mardec.remove_coord('time')
            #
            # Loop over remaining years and add contributions
            for yearc in range(self.year1+1,self.year2+1):
                # Do not use 1978 for olrinterp as missing data from 16 Mar 1978
                # to 31 Dec 1978.  Just use data from 1979 onwards to calculate
                # annual cycle
                if self.data_source=='olrinterp' and yearc==1978:
                    raise UserWarning('Cannot use 1978 for olrinterp because of missing data.')
                # Extract 1 Jan to 28 Feb of current year
                time1=datetime.datetime(yearc,1,1)
                time2=datetime.datetime(yearc,2,28)
                print(time1,time2)
                time_constraint=iris.Constraint(time=lambda cell: time1<=cell<=time2)
                xx1=self.data_in.extract(time_constraint)
                janfebc=xx1.concatenate_cube()
                janfeb_counterc,janfebc=create_counter_from_mask(janfebc,verbose=self.verbose)
                janfebc.remove_coord('time')
                # Extract 1 Mar to 31 Dec of current year
                time1=datetime.datetime(yearc,3,1)
                time2=datetime.datetime(yearc,12,31)
                print(time1,time2)
                time_constraint=iris.Constraint(time=lambda cell: time1<=cell<=time2)
                xx1=self.data_in.extract(time_constraint)
                mardecc=xx1.concatenate_cube()
                mardec_counterc,mardecc=create_counter_from_mask(mardecc,verbose=self.verbose)
                mardecc.remove_coord('time')
                # Add contributions for current year
                janfeb_counter+=janfeb_counterc
                mardec_counter+=mardec_counterc
                janfeb=iris.analysis.maths.add(janfeb,janfebc)
                mardec=iris.analysis.maths.add(mardec,mardecc)
            # Save units and copy them back later
            units=janfeb.units
            # Create and add a time coordinate for year yearac, and attributes
            time_first=datetime.datetime(year=self.yearac,month=1,day=1)
            time_units='days since '+str(time_first)
            # janfeb
            time_val=np.arange(59)
            tcoord=iris.coords.DimCoord(time_val,standard_name='time',units=time_units)
            janfeb.add_dim_coord(tcoord,0)
            janfeb.rename(janfebc.name())
            janfeb.var_name=janfebc.var_name
            janfeb.attributes=janfebc.attributes
            # mardec
            time_val=np.arange(59,365)
            tcoord=iris.coords.DimCoord(time_val,standard_name='time',units=time_units)
            mardec.add_dim_coord(tcoord,0)
            mardec.rename(mardecc.name())
            mardec.var_name=mardecc.var_name
            mardec.attributes=mardecc.attributes
            # Recombine janfeb and mardec
            xx1=iris.cube.CubeList([janfeb,mardec])        
            xx1=xx1.concatenate_cube()
            counter=np.concatenate((janfeb_counter,mardec_counter))
        elif self.calendar=='360_day':
            # First year
            yearc=self.year1
            # Extract 1 Jan to 30 Dec of current year
            time1=cftime.Datetime360Day(yearc,1,1)
            time2=cftime.Datetime360Day(yearc,12,30)
            print(time1,time2)
            time_constraint=iris.Constraint(time=lambda cell: time1<=cell.point<=time2)
            xx1=self.data_in.extract(time_constraint)
            allyear=xx1.concatenate_cube()
            allyear_counter,allyear=create_counter_from_mask(allyear,verbose=self.verbose)
            allyear.remove_coord('time')
            # Loop over remaining years and add contributions
            for yearc in range(self.year1+1,self.year2+1):
                # Extract 1 Jan to 30 Dec of current year
                time1=cftime.Datetime360Day(yearc,1,1)
                time2=cftime.Datetime360Day(yearc,12,30)
                print(time1,time2)
                time_constraint=iris.Constraint(time=lambda cell: time1<=cell.point<=time2)
                xx1=self.data_in.extract(time_constraint)
                allyearc=xx1.concatenate_cube()
                allyear_counterc,allyearc=create_counter_from_mask(allyearc,verbose=self.verbose)
                allyearc.remove_coord('time')
                # Add contributions for current year
                allyear_counter+=allyear_counterc
                allyear=iris.analysis.maths.add(allyear,allyearc)
            # Save units and copy them back later
            units=allyear.units
            # Create and add a time coordinate for year yearac, and attributes
            time_first=cftime.Datetime360Day(year=self.yearac,month=1,day=1)
            time_units=cf_units.Unit('days since '+str(time_first),calendar=self.calendar)
            # allyear
            time_val=np.arange(360)
            tcoord=iris.coords.DimCoord(time_val,standard_name='time',units=time_units)
            allyear.add_dim_coord(tcoord,0)
            allyear.rename(allyearc.name())
            allyear.var_name=allyearc.var_name
            allyear.attributes=allyearc.attributes
            xx1=allyear
            counter=allyear_counter
        else:
            raise UserWarning('calendar not recognised.')
        counter_min=counter.min()
        counter_max=counter.max()
        counter_mean=counter.mean()
        print('counter_min',counter_min)
        print('kmin',self.kmin)
        print('counter_max',counter_max)
        print('counter_mean',counter_mean)

        # If there are grid points in the counter arrays that have a count
        # less than kmin, create a masked counter array
        if counter_min<self.kmin:
            print('Some grid points have < kmin values. Adding a mask to counter.')
            lmask=True
            mask=np.less(counter,self.kmin)
            counter=np.ma.masked_array(counter,mask=mask)
            # Check there are not too many missing grid points
            npts=np.prod(mask.shape)
            npts_missing=mask.sum()
            frac_missing=npts_missing/npts
            print('npts_missing,npts',npts_missing,npts)
            print('frac_missing,frac_crit',frac_missing,self.frac_crit)
            if (frac_missing>self.frac_crit):
                raise UserWarning('Warning: Too many missing points. Rethink algorithm.')
        else:
            print('All grid points have > kmin values. No mask needed for counter.')
            lmask=False
        xx2=iris.analysis.maths.divide(xx1,counter)
        xx1=create_cube(conv_float32(xx2.data),xx1)
        # Replace original units
        xx1.units=units
        # Set to mean where less than kmin data points went in to calculation
        mean=xx1.data.mean()
        print('mean',mean)
        if counter_min<self.kmin:
            xx2=np.where(mask,mean,xx1.data)
            xx1=create_cube(conv_float32(xx2),xx1)
        # Add a cell method to describe the action of creating annual cycle
        comments='raw annual cycle '+str(self.year1)+'-'+str(self.year2)+': kmin='+str(self.kmin)+': frac_crit='+str(self.frac_crit)
        cm=iris.coords.CellMethod('mean','time',comments=comments)
        xx1.add_cell_method(cm)
        # Set data_anncycle_raw attribute
        self.data_anncycle_raw=xx1
        # Save raw annual cycle
        iris.save(self.data_anncycle_raw,self.file_anncycle_raw)
        if self.archive:
            archive_file(self,self.file_anncycle_raw)

    def f_read_anncycle_raw(self):
        """Read previously created raw annual cycle.

        Create data_anncycle_raw attribute.
        
        """
        self.data_anncycle_raw=iris.load_cube(self.file_anncycle_raw,self.name)

    def f_anncycle_smooth(self):
        """Calculate smoothed annual cycle.

        This is created from the mean and first self.nharm annual
        harmonics of the raw annual cycle.

        f(t_i)= \overline f + \Sum_{k=1}^\{nharm} A_k \cos \omega_k t
                                                + B_k \sin \omega_k t

        \omega_k = 2\pi k / T

        A_k = 2/N \Sum_{k=1}^\{N} f(t_i) \cos \omega_k t_i

        B_k = 2/N \Sum_{k=1}^\{N} f(t_i) \sin \omega_k t_i

        Create data_mean, data_anncycle_smooth attributes.
        Create data_anncycle_smooth_leap attribute if calendar is gregorian.
        
        """
        if self.frequency!='d':
            raise ToDoError('Annual cycle presently only coded up for daily data.')
        # Calculate annual mean
        print('data_anncycle_raw.shape: {0!s}'.format(self.data_anncycle_raw.shape))
        self.data_mean=self.data_anncycle_raw.collapsed('time',iris.analysis.MEAN)
        print('data_mean.shape: {0!s}'.format(self.data_mean.shape))
        # Calculate cosine and sine harmonics of annual cycle
        # Create a (ntime,1) array of itime
        ntime=self.data_anncycle_raw.coord('time').shape[0]
        itime=np.array(np.arange(ntime))
        itime=itime.reshape(itime.shape+(1,))
        print('itime.shape: {0!s}'.format(itime.shape))
        # Create a (1,nharm) array of iharm
        iharm=np.array(np.arange(1,self.nharm+1))
        iharm=iharm.reshape(iharm.shape+(1,))
        iharm=iharm.transpose()
        print('iharm: {0!s}'.format(iharm))
        print('iharm.shape: {0!s}'.format(iharm.shape))
        # Create (ntime,nharm) arrays of cosine and sine waves
        argument=(2*np.pi/float(ntime))*np.dot(itime,iharm)
        cosine=np.cos(argument)
        sine=np.sin(argument)
        print('argument: {0!s}'.format(argument))
        print('cosine: {0!s}'.format(cosine))
        print('cosine.shape: {0!s}'.format(cosine.shape))
        # Duplicate data array with an nharm-length axis and reshape to
        # (nlat,nlon,...,ntime,nharm) in xx1
        # NB the nlat,nlon,... could be just nlat,nlon, or eg nlat,lon,nlev
        xx1=self.data_anncycle_raw
        print('Original order. xx1.shape: {0!s}'.format(xx1.shape))
        # Find index of time coordinate
        dim_coord_names=[xx.name() for xx in xx1.dim_coords]
        time_index=dim_coord_names.index('time')
        print('time_index',time_index)
        # Reorder such that time is final index
        ndim=len(xx1.shape)
        print('ndim',ndim)
        xx1=np.moveaxis(xx1.data,time_index,ndim-1)
        print('Order should be (...,t). xx1.shape: {0!s}'.format(xx1.shape))
        # Duplicate data and create extra dimension of length nharm
        xx1=xx1.reshape(xx1.shape+(1,))
        xx1=np.dot(xx1,np.ones((1,self.nharm)))
        print('Order should be (...,t,nharm). xx1.shape: {0!s}'.format(xx1.shape))
        xx2=xx1
        xx2cos=xx2*cosine
        xx2sin=xx2*sine
        print('Time should be index -2. xx2cos.shape: {0!s}'.format(xx2cos.shape))
        # Sum over time axis (index is -2)
        xx2cos=xx2cos.sum(axis=-2)
        xx2sin=xx2sin.sum(axis=-2)
        # Multiply by 2/ntime to get (nlat,nlon,...,nharm) array of a and b
        # A_k and B_k coefficients
        a_coeff=2.0*xx2cos/float(ntime)
        b_coeff=2.0*xx2sin/float(ntime)
        print('a_coeff.shape: {0!s}'.format(a_coeff.shape))
        # Matrix multiply the a,b coefficients by the (nharm,ntime) arrays of
        # cosine and sine harmonic time series to get the contributions to the
        # smoothed annual cycle in a (nlat,nlon,...,ntime) array
        xx2cos=np.dot(a_coeff,cosine.transpose())
        xx2sin=np.dot(b_coeff,sine.transpose())
        print('Time should be last dimension. xx2cos.shape: {0!s}'.format(xx2cos.shape))
        # Reshape so it can be added to mean by broadcasting
        xx3cos=np.moveaxis(xx2cos,-1,0)
        xx3sin=np.moveaxis(xx2sin,-1,0)
        print('Time should be first dimension. xx3cos.shape: {0!s}'.format(xx3cos.shape))
        # Add cosine and sine contributions and time mean to get smoothed
        # annual cycle
        # NB Only add time mean if not detrending, otherwise time mean would be removed twice!
        if self.detrend:
            xx4=xx3cos+xx3sin
        else:
            xx4=self.data_mean.data+xx3cos+xx3sin
        # Create a cube from this numpy array
        xx4=create_cube(conv_float32(xx4),self.data_anncycle_raw)
        # Add a cell method to describe the smoothed annual cycle
        cm=iris.coords.CellMethod('mean','time',comments='smoothed annual cycle: mean + '+str(self.nharm)+' harmonics')
        xx4.add_cell_method(cm)
        # Set data_anncycle_smooth attribute
        self.data_anncycle_smooth=xx4
        print('data_anncycle_smooth',self.data_anncycle_smooth)
        # Save smoothed annual cycle
        iris.save(self.data_anncycle_smooth,self.file_anncycle_smooth)
        if self.archive:
            archive_file(self,self.file_anncycle_smooth)
        #
        if self.calendar=='gregorian':
            print('# Create alternate smoothed annual cycle for leap year, with 29 Feb')
            print('# represented by a copy of 28 Feb.')
            # Use year self.yearacleap for this, as this is a leap year.
            time_units_leap='days since '+str(self.yearacleap).zfill(4)+'-01-01 00:00:0.0'
            # 1 Jan to 28 Feb of smoothed annual cycle
            time1=datetime.datetime(self.yearac,1,1)
            time2=datetime.datetime(self.yearac,2,28)
            time_constraint=iris.Constraint(time=lambda cell: time1<=cell<=time2)
            janfeb=xx4.extract(time_constraint)
            time_val=janfeb.coord('time').points
            tcoord2=iris.coords.DimCoord(time_val,standard_name='time',var_name='time',units=time_units_leap)
            janfeb.remove_coord('time')
            janfeb.add_dim_coord(tcoord2,0)
            # 29 Feb only of smoothed annual cycle. Copy of 28 Feb data.
            time1=datetime.datetime(self.yearac,2,28)
            time_constraint=iris.Constraint(time=time1)
            feb29=xx4.extract(time_constraint)
            feb29=iris.util.new_axis(feb29,'time') # Promote singleton aux coord to dim coord
            time_val=np.array([59,]) # Julian day number for 29 Feb
            tcoord2=iris.coords.DimCoord(time_val,standard_name='time',var_name='time',units=time_units_leap)
            feb29.remove_coord('time')
            feb29.add_dim_coord(tcoord2,0)
            # 1 Mar to 31 Dec of smoothed annual cycle
            time1=datetime.datetime(self.yearac,3,1)
            time2=datetime.datetime(self.yearac,12,31)
            time_constraint=iris.Constraint(time=lambda cell: time1<=cell<=time2)
            mardec=xx4.extract(time_constraint)
            time_val=mardec.coord('time').points+1
            tcoord2=iris.coords.DimCoord(time_val,standard_name='time',var_name='time',units=time_units_leap)
            mardec.remove_coord('time')
            mardec.add_dim_coord(tcoord2,0)
            # Create a cubelist then a single cube
            xx1=iris.cube.CubeList([janfeb,feb29,mardec])
            xx1=xx1.concatenate_cube()
            # Set data_anncycle_smooth_leap attribute
            self.data_anncycle_smooth_leap=xx1
            print('data_anncycle_smooth_leap',self.data_anncycle_smooth_leap)
            # Save smoothed annual cycle for leap year
            iris.save(self.data_anncycle_smooth_leap,self.file_anncycle_smooth_leap)
            if self.archive:
                archive_file(self,self.file_anncycle_smooth_leap)

    def f_read_anncycle_smooth(self):
        """Read previously calculated smoothed annual cycle.

        Create data_anncycle_smooth and data_anncycle_smooth_leap
        attributes.
        
        """
        self.data_anncycle_smooth=iris.load_cube(self.file_anncycle_smooth,self.name)
        if self.calendar=='gregorian':
            self.data_anncycle_smooth_leap=iris.load_cube(self.file_anncycle_smooth_leap,self.name)

    def f_expand_anncycle_smooth(self):
        """Expand smoothed annual cycle of daily data for higher frequency input data.

        For e.g., Gregorian non-leap year, the daily data smoothed
        annual cycle has 365 days. If the input data is '6h', then the
        expanded annual cycle will have 365*4=1460 days, with the
        daily values repeated as necessary. Here, the 1 Jan data will
        be copied 4 times, for 0000, 0600, 1200, 1800 on 1
        Jan. Similarly for 2 Jan, etc.

        Create data_expanded_anncycle_smooth attribute, and
        data_expanded_anncycle_smooth_leap attribute if calendar is
        Gregorian.

        """
        if not self.anncycle_source:
            raise UserWarning('Only call this function if using higher than daily frequency input data.')
        npd=find_npd(self.source)
        fracday=1/npd
        print('npd,fracday: {0!s}, {1!s}'.format(npd,fracday))
        def expanded_anncycle(anncycle_d,npd):
            tcoord_d=anncycle_d.coord('time')
            tcoord_hf_vals=[]
            for dayc in tcoord_d.points:
                for ii in range(npd):
                    xx=dayc+ii*fracday
                    tcoord_hf_vals.append(xx)
            tcoord_hf=iris.coords.DimCoord(tcoord_hf_vals,standard_name='time',units=tcoord_d.units)
            shape_d=anncycle_d.shape
            shape_hf=(shape_d[0]*npd,)+shape_d[1:]
            print('shape_d,shape_hf: {0!s}, {1!s}'.format(shape_d,shape_hf))
            x1=1e20*np.ones(shape_hf) # create empty array then overwrite
            ntime_d=tcoord_d.shape[0]
            for itime in range(ntime_d):
                xx=anncycle_d.data[itime]
                for ii in range(npd):
                    kk=itime*npd+ii
                    print(itime,ii,kk)
                    x1[kk]=xx
            x2=create_cube(conv_float32(x1),anncycle_d,new_axis=tcoord_hf)
            x2.attributes['frequency']=self.frequency
            cm=iris.coords.CellMethod('mean','time',comments='Smoothed annual cycle: daily expanded to '+self.frequency)
            x2.add_cell_method(cm)
            return x2
        # Standard year
        self.data_expanded_anncycle_smooth=expanded_anncycle(self.data_anncycle_smooth,npd)
        iris.save(self.data_expanded_anncycle_smooth,self.file_expanded_anncycle_smooth)
        if self.archive:
            archive_file(self,self.file_expanded_anncycle_smooth)
        # Leap year
        if self.calendar=='gregorian':
            self.data_expanded_anncycle_smooth_leap=expanded_anncycle(self.data_anncycle_smooth_leap,npd)
            iris.save(self.data_expanded_anncycle_smooth_leap,self.file_expanded_anncycle_smooth_leap)
            if self.archive:
                archive_file(self,self.file_expanded_anncycle_smooth_leap)
        # Overwrite data_anncycle_smooth and data_anncycle_smooth_leap with expanded versions
        self.data_anncycle_smooth=self.data_expanded_anncycle_smooth
        if self.calendar=='gregorian':
            self.data_anncycle_smooth_leap=self.data_expanded_anncycle_smooth_leap

    def f_read_expanded_anncycle_smooth(self):
        """Read previously calculated expanded smoothed annual cycle.

        Create data_anncycle_smooth and data_anncycle_smooth_leap
        attributes.
        
        """
        self.data_anncycle_smooth=iris.load_cube(self.file_expanded_anncycle_smooth,self.name)
        if self.calendar=='gregorian':
            self.data_anncycle_smooth_leap=iris.load_cube(self.file_expanded_anncycle_smooth_leap,self.name)

    def f_subtract_anncycle(self):
        """Subtract (trend and) smoothed annual cycle from input data to create anomaly data.

        Read in input data and process and write output data (anomaly
        with (trend and) smoothed annual cycle subtracted) in chunks of
        outfile_frequency, e.g., 'year'.

        Makes allowance for input data not having to start on 1 Jan or
        end on 31 Dec.

        Create data_anncycle_rm attribute.
        
        """
        # Set initial value of current year
        yearc=self.time1.year
        # Set final value of current year
        year_end=self.time2.year
        # Loop over years
        while yearc<=year_end:
            # If current year is leap year and gregorian calendar, use
            # self.data_anncycle_smooth_leap, otherwise use
            # self.data_anncycle_smooth.
            if divmod(yearc,4)[1]==0 and self.calendar=='gregorian':
                leap=True
                smooth_year_number=self.yearacleap
                anncycle_smooth=self.data_anncycle_smooth_leap
            else:
                leap=False
                smooth_year_number=self.yearac
                anncycle_smooth=self.data_anncycle_smooth
            print('yearc,leap: {0!s}, {1!s}'.format(yearc,leap))
            # Set start and end times for input and anncycle data
            if self.calendar=='gregorian':
                time_beg_in=datetime.datetime(yearc,1,1,0,0)
                if time_beg_in<self.time1:
                    time_beg_in=self.time1
                time_end_in=datetime.datetime(yearc,12,31,23,59)
                if time_end_in>self.time2:
                    time_end_in=self.time2
                time_beg_anncycle=datetime.datetime(smooth_year_number,time_beg_in.month,time_beg_in.day,0,0)
                time_end_anncycle=datetime.datetime(smooth_year_number,time_end_in.month,time_end_in.day,23,59)
            elif self.calendar=='360_day':
                time_beg_in=cftime.Datetime360Day(yearc,1,1,0,0)
                if time_beg_in<self.time1:
                    time_beg_in=self.time1
                time_end_in=cftime.Datetime360Day(yearc,12,30,23,59)
                if time_end_in>self.time2:
                    time_end_in=self.time2
                time_beg_anncycle=cftime.Datetime360Day(smooth_year_number,time_beg_in.month,time_beg_in.day,0,0)
                time_end_anncycle=cftime.Datetime360Day(smooth_year_number,time_end_in.month,time_end_in.day,23,59)
            else:
                raise UserWarning('Invalid calendar.')
            print(time_beg_in,time_end_in,time_beg_anncycle,time_end_anncycle)
            # Extract input data for this (potentially partial) year
            time_constraint=set_time_constraint(time_beg_in,time_end_in,calendar=self.calendar,verbose=self.verbose)
            x1=self.data_in.extract(time_constraint)
            data_in=x1.concatenate_cube()
            # Extract anncycle data for this (potentially partial) year
            time_constraint=set_time_constraint(time_beg_anncycle,time_end_anncycle,calendar=self.calendar,verbose=self.verbose)
            data_anncycle=anncycle_smooth.extract(time_constraint)
            # Remove time coordinates from input and anncycle cubes so they
            # can be subtracted
            tcoord=data_in.coord('time')
            data_in.remove_coord('time')
            data_anncycle.remove_coord('time')
            if self.detrend:
                # Calculate trend terms for this time chunk and subtract
                yy=np.dot(tcoord.points.reshape(tcoord.points.shape+(1,)),self.mm.data.reshape((1,self.ngrid)))
                print('yy.shape,tcoord.shape,mm.shape: {0!s}, {1!s}, {2!s}'.format(yy.shape,tcoord.points.shape,self.mm.shape))
                yy=yy.reshape(data_in.shape)
                yy=yy+self.cc.data
                print('yy.shape: {0!s}'.format(yy.shape))
                yy=create_cube(yy,data_in,new_var_name='dummy')
                xx1=iris.analysis.maths.subtract(data_in,yy)
            else:
                xx1=data_in
            # Subtract smoothed annual cycle for current period to create anomaly
            data_anom=iris.analysis.maths.subtract(xx1,data_anncycle)
            # Add back time coordinate and update metadata
            data_anom.add_dim_coord(tcoord,0)
            data_anom.rename(data_in.name())
            data_anom.var_name=data_in.var_name
            data_anom.attributes=data_in.attributes
            # Add a cell method to describe the smoothed annual cycle
            # There is a small list of allowed method names (1st argument)
            # 'point' seems the most appropriate, actions on each point?!
            if self.detrend:
                str1='anomaly: subtracted linear trend and smoothed annual cycle: '
            else:
                str1='anomaly: subtracted smoothed annual cycle: '
            cm=iris.coords.CellMethod('point','time',comments=str1+str(self.year1)+'-'+str(self.year2)+': mean + '+str(self.nharm)+' harmonics')
            data_anom.add_cell_method(cm)
            # Save anomaly data (anncycle subtracted)
            self.year=yearc
            if self.outfile_frequency=='year':
                fileout=replace_wildcard_with_time(self,self.file_anncycle_rm)
                print('fileout: {0!s}'.format(fileout))
                iris.save(data_anom,fileout)
                if self.archive:
                    archive_file(self,fileout)
            elif self.outfile_frequency=='month':
                month1=time_beg_in.month
                month2=time_end_in.month
                for monthc in range(time_beg_in.month,time_end_in.month+1):
                    self.month=monthc
                    fileout=replace_wildcard_with_time(self,self.file_anncycle_rm)
                    print('fileout: {0!s}'.format(fileout))
                    if self.calendar=='gregorian':
                        time1=datetime.datetime(yearc,monthc,1)
                        if monthc!=12:
                            time2=datetime.datetime(yearc,monthc+1,1)
                        else:
                            time2=datetime.datetime(yearc+1,1,1)
                        time2-=datetime.timedelta(days=1)
                    elif self.calendar=='360_day':
                        time1=cftime.Datetime360Day(yearc,monthc,1)
                        time2=cftime.Datetime360Day(yearc,monthc,30)
                    print('time1,time2: {0!s}, {1!s}'.format(time1,time2))
                    time_constraint=set_time_constraint(time1,time2,calendar=self.calendar,verbose=self.verbose)
                    x1=data_anom.extract(time_constraint)
                    iris.save(x1,fileout)
                    if self.archive:
                        archive_file(self,fileout)
            else:
                raise ToDoError('Need code for other outfile_frequency values.')
            # Increment current year
            yearc+=1
        # Set data_anncycle_rm attribute
        self.data_anncycle_rm=iris.load(self.file_anncycle_rm,self.name)

    def f_read_subtract_anncycle(self):
        """Read previously calculatedanomaly data (annual cycle subtracted).

        Create data_anncycle_rm attribute.
        
        """
        self.data_anncycle_rm=iris.load(self.file_anncycle_rm,self.name)
        
#==========================================================================

class GliderMission(object):
    
    """Analyse a glider mission

    Called from, e.g., glider_interp_lon.py

    Selected attributes:

    self.mission : integer id for mission

    self.gliderids : list of integer glider ids for mission

    self.time1 : datetime.datetime object for 00 UTC on day of
    deployment of first glider in mission

    self.time2 : datetime.datetime object for 00 UTC on day after
    recovery of last glider in mission

    """

    def __init__(self,**descriptor):
        """Initialise from descriptor dictionary.

        Compulsory keywords: 'verbose','source_wildcard','var_name',
        'basedir','mission'.
        """
        self.__dict__.update(descriptor)
        self.descriptor=descriptor
        self.name=var_name2long_name[self.var_name]
        if self.mission==31:
            self.gliderids=[579,534,532,620,613]
            self.time1=datetime.datetime(2016,6,30)
            self.time2=datetime.datetime(2016,7,21)
        else:
            raise UserWarning('Invalid mission')
        # self.gliders is a dictionary of glider objects
        self.gliders={}
        for gliderid in self.gliderids:
            self.gliders[gliderid]=Glider(gliderid,**descriptor)
        self.data_oi_pad_all={}
        self.data_oi_interp_lon={}
        if self.verbose:
            print(self)        

    def __repr__(self):
        return 'GliderMission({0.descriptor!r},verbose={0.verbose!r})'.format(self)

    def __str__(self):
        if self.verbose==2:
            ss=h1a+'GliderMission instance \n'+\
                'mission: {0.mission!s} \n'+\
                'gliderids: {0.gliderids!s} \n'+\
                'source_wildcard: {0.source_wildcard!s} \n'+\
                'time1: {0.time1!s} \n'+\
                'time2: {0.time2!s} \n'+h1b
            return(ss.format(self))
        else:
            return 'GliderMission instance'

    def f_interp_oi_lon(self,lon1,lon2,delta_lon):
        """Interpolate padded OI data from individual gliders in longitude."""

        x1a=iris.cube.CubeList([]) # var_name
        x2a=iris.cube.CubeList([]) # longitude
        kount=0
        for gliderid in self.gliderids:
            print('gliderid: {0!s}'.format(gliderid))
            # Read and pad OI var_name
            self.gliders[gliderid].f_read_oi(self.var_name,verbose=False)
            self.gliders[gliderid].f_oi_pad(self.var_name,self.time1,self.time2,verbose=False)
            # Read and pad OI longitude
            self.gliders[gliderid].f_read_oi('lon',verbose=False)
            self.gliders[gliderid].f_oi_pad('lon',self.time1,self.time2,verbose=False)
            # Combine padded OI data from all gliders into single iris cube
            glider_coord=iris.coords.DimCoord(kount,var_name='glider')
            # Variable
            x1b=self.gliders[gliderid].data_oi_pad[self.var_name]
            x1b.add_aux_coord(glider_coord)
            x1b.data.fill_value=1e20 # Reset or merge_cube will fail
            x1a.append(x1b)
            # Longitude
            x2b=self.gliders[gliderid].data_oi_pad['lon']
            x2b.add_aux_coord(glider_coord)
            x2b.data.fill_value=1e20 # Reset or merge_cube will fail
            x2a.append(x2b)
            kount+=1
        # Merge cube list of individual gliders to cube of all gliders
        x1c=x1a.merge_cube()
        x2c=x2a.merge_cube()
        # Create data_oi_pad_all attribute dictionary entries
        self.data_oi_pad_all[self.var_name]=x1c
        self.data_oi_pad_all['lon']=x2c
        # Create longitude axis to interpolate onto
        delta=1e-6
        xnew=np.arange(lon1,lon2+delta,delta_lon)
        print('xnew: {0!s}'.format(xnew))
        nlon=len(xnew)
        lon_coord=iris.coords.DimCoord(xnew,standard_name='longitude',var_name='longitude',units='degree_east')
        print('lon_coord: {0!s}'.format(lon_coord))
        # Interpolation
        xx=self.data_oi_pad_all['lon'].data
        yy=self.data_oi_pad_all[self.var_name].data
        # As the longitudes (x values) of each glider are different at each
        # time, depth etc, have to do a separate interpolation over longitude
        # at each time, depth etc.  Ugly nested loops
        if len(yy.shape)!=3:
            raise UserWarning('Need some more code if not 3-d data')
        # Create empty new 3-d array to fill with interpolated values
        shape=(nlon,)+yy.shape[1:]
        ndim1=shape[1]
        ndim2=shape[2]
        print('nlon,ndim1,ndim2: {0!s},{1!s},{2!s}'.format(nlon,ndim1,ndim2))
        x3=np.zeros(shape)
        # Loop over other dimensions and fill x3 with interpolated values
        missing_value=1e20
        nglider=len(self.gliderids)
        for idim1 in range(ndim1):
            print('idim1: {0!s}'.format(idim1))
            for idim2 in range(ndim2):
                xx1=xx[:,idim1,idim2]
                yy1=yy[:,idim1,idim2]
                # Get rid of masked data
                xxc=[xx1.data[ii] for ii in range(nglider) if not xx1.mask[ii]]
                yyc=[yy1.data[ii] for ii in range(nglider) if not xx1.mask[ii]]
                n_data_points=len(xxc)
                #print('idim1,idim2,n_data_points: {0!s},{1!s},{2!s}'.format(idim1,idim2,n_data_points))
                #print('xxc: {0!s}'.format(xxc))
                #print('yyc: {0!s}'.format(yyc))
                if n_data_points>1:
                    # Interpolate
                    ynew=np.interp(xnew,xxc,yyc,left=missing_value,right=missing_value)
                elif n_data_points==1:
                    # Only 1 data point.  
                    # Take the data at the 1 point and
                    # extend it over neighbouring points.
                    xspread=3*delta_lon
                    ynew=np.where(abs(xxc[0]-xnew)<xspread,yyc[0],missing_value)
                else:
                    # Zero data points.
                    # Fill ynew with missing values
                    ynew=missing_value*np.ones((nlon,))
                #print('ynew: {0!s}'.format(ynew))
                x3[:,idim1,idim2]=ynew
        # Create iris cube
        oldcube=self.data_oi_pad_all[self.var_name]
        # Create a list of two-lists, each of form [dim_coord,index]
        dim_coords=[[lon_coord,0]]
        kdim=1
        for xx in oldcube.dim_coords[1:]:
            dim_coords.append([xx,kdim])
            kdim+=1
        x4=iris.cube.Cube(conv_float32(x3),units=oldcube.units,attributes=oldcube.attributes,cell_methods=oldcube.cell_methods,dim_coords_and_dims=dim_coords)
        x4.rename(oldcube.name())
        x4.var_name=oldcube.var_name
        # Mask missing values
        x4.data=np.ma.masked_greater(x4.data,missing_value/10)
        # Create data_oi_interp_lon attribute dictionary entry
        self.data_oi_interp_lon[self.var_name]=x4
        # Save interpolated data
        source_out=self.source_wildcard.replace('???','all')
        #
        ########################################################
        # Need to sort out code with '2016' below to handle year properly
        # in gliderMission and Glider classes
        ########################################################
        #
        fileout=os.path.join(self.basedir,source_out,'std',self.var_name+'_all_2016.nc')
        iris.save(self.data_oi_interp_lon[self.var_name],fileout)
        if self.archive:
            archive_file(self,fileout)
        
        
#==========================================================================

class Glider(object):

    """Glider object."""

    def __init__(self,gliderid,**descriptor):
        """Initialise from descriptor dictionary.

        Compulsory keywords: 'verbose','source_wildcard',
        'basedir','mission'.
        """
        self.__dict__.update(descriptor)
        self.descriptor=descriptor
        self.gliderid=gliderid
        self.data_oi={}
        self.data_oi_pad={}
        self.verbose=verbose

    def __repr__(self):
        return 'Glider ({0.gliderid!r},verbose={0.verbose!r})'.format(self)

    def __str__(self):
        if self.verbose==2:
            ss=h1a+'Glider instance \n'+\
                'gliderid: {0.gliderid!s} \n'+\
                'mission: {0.mission!s} \n'+\
                'source_wildcard: {0.source_wildcard!s} \n'+h1b
            return(ss.format(self))
        else:
            return 'Glider instance'

    def f_read_oi(self,var_name,verbose=False):
        """Read OI field for individual glider.

        Create attribute data_oi.
        """
        name=var_name2long_name[var_name]
        self.source=self.source_wildcard.replace('???',str(self.gliderid))
        file1=os.path.join(self.basedir,self.source,'std',var_name+'_all_????.nc')
        x1=iris.load(file1,name)
        x2=x1.concatenate_cube()
        # Set data_oi attribute dictionary entry
        self.data_oi[var_name]=x2
        if verbose:
            print('source {0.source!s}'.format(self))

    def f_oi_pad(self,var_name,time1,time2,verbose=False):
        """Pad OI data with missing values back to time1 and forward to time2.

        Create attribute data_oi_pad.
        """
        name=var_name2long_name[var_name]
        source_info(self)
        # Get first and last time of data: time1a,time2a
        tcoord=self.data_oi[var_name].coord('time')
        time_units=tcoord.units
        time1a=tcoord.units.num2date(tcoord.points[0])
        time2a=tcoord.units.num2date(tcoord.points[-1])
        time1a_val=time_units.date2num(time1a)
        time2a_val=time_units.date2num(time2a)
        if verbose:
            print('time1a,time2a: {0!s}, {1!s}'.format(time1a,time2a))
            print('time1a_val,time2a_val: {0!s}, {1!s}'.format(time1a_val,time2a_val))
        # Process times to pad to: time1,time2
        time1_val=time_units.date2num(time1)
        time2_val=time_units.date2num(time2)
        if verbose:
            print('time1,time2: {0!s}, {1!s}'.format(time1,time2))
            print('time1_val,time2_val: {0!s}, {1!s}'.format(time1_val,time2_val))
        if time1>time1a or time2<time2a:
            raise UserWarning('time1 or time2 not set correctly')
        # Find index of time coordinate
        missing_value=1e20
        kount=0
        for dimc in self.data_oi[var_name].dim_coords:
            if dimc.standard_name=='time':
                tcoord_index=kount
                kount+=1
        if verbose:
            print('tcoord_index: {0!s}'.format(tcoord_index))
        #
        # Create padding array with missing values between time1 and time1a
        ntime=int((time1a-time1)/self.timedelta)
        shape=self.data_oi[var_name].shape
        shape1=shape[:tcoord_index]+(ntime,)+shape[tcoord_index+1:]
        if verbose:
            print('ntime: {0!s}'.format(ntime))
            print('shape: {0!s}'.format(shape))
            print('shape1: {0!s}'.format(shape1))
        x1=missing_value*np.ones(shape1)
        tcoord1_vals=[time_units.date2num(time1+ii*self.timedelta) for ii in range(ntime)]
        tcoord1_vals=conv_float32(np.array(tcoord1_vals))
        tcoord1=iris.coords.DimCoord(tcoord1_vals,standard_name='time',var_name='time',units=time_units)        
        if verbose:
            print('tcoord1_vals: {0!s}'.format(tcoord1_vals))
            print('tcoord1: {0!s}'.format(tcoord1))
        pad_before=create_cube(conv_float32(x1),self.data_oi[var_name],new_axis=tcoord1)
        #
        # Create padding array with missing values between time2a and time2'
        ntime=int((time2-time2a)/self.timedelta)
        shape=self.data_oi[var_name].shape
        shape1=shape[:tcoord_index]+(ntime,)+shape[tcoord_index+1:]
        if verbose:
            print('ntime: {0!s}'.format(ntime))
            print('shape: {0!s}'.format(shape))
            print('shape1: {0!s}'.format(shape1))
        x1=missing_value*np.ones(shape1)
        tcoord1_vals=[time_units.date2num(time2a+(ii+1)*self.timedelta) for ii in range(ntime)]
        tcoord1_vals=conv_float32(np.array(tcoord1_vals))
        tcoord1=iris.coords.DimCoord(tcoord1_vals,standard_name='time',var_name='time',units=time_units)        
        if verbose:
            print('tcoord1_vals: {0!s}'.format(tcoord1_vals))
            print('tcoord1: {0!s}'.format(tcoord1))
        pad_after=create_cube(conv_float32(x1),self.data_oi[var_name],new_axis=tcoord1)
        #
        # Create cube list of padding before, data, padding after
        x4=iris.cube.CubeList([pad_before,self.data_oi[var_name],pad_after])
        x5=x4.concatenate_cube()
        # Mask missing values
        x5.data=np.ma.masked_greater(x5.data,missing_value/10)
        # Set data_oi_pad attribute dictionary entry
        self.data_oi_pad[var_name]=x5
        
#==========================================================================

class CubeDiagnostics(object):

    """CubeDiagnostics object.

    An object to calculate physical diagnostics of data, e.g.,
    calculate the mixed layer depth using a particular method, from a
    cube of conservative temperature.

    NB There are no generic diagnostics here, such as calculation of
    the x-derivative (which could be applied to a cube of any
    spatially varying quantity).  Instead, reserve this class to
    contain methods that only work on specific physical quantities,
    such as the above example.

    Other examples might be calculating density from pressure and
    temperature using the ideal gas law.

    An instance of the class can have attributes such as self.tsc
    (conservative_temperature), self.rho (air_density), etc.  These
    will be iris cubes.  Each iris cube should have the same dimensions.
    
    After creating an instance of the class, call method f_read_data()
    to lazy read whatever variables (e.g., 'tsc', or 'rho' and 'ta')
    will be needed later.

    Then loop over relevant time blocks in main programme, and call
    desired method from this class, e.g., f_mld().  These methods,
    e.g., f_mld() then read and set attributes named after the
    variable name, e.g., self.tsc is an iris cube of tsc, and
    calculate diagnostics, e.g., self.mld is an iris cube of mixed
    layer depth.

    It is anticipated that this class will grow over time and will
    eventually be very large, adding more methods as needed.

    """

    def __init__(self,**descriptor):
        """Initialise from descriptor dictionary.

        Compulsory keywords: 'verbose','source','level'
        'basedir'.
        """
        self.__dict__.update(descriptor)
        self.descriptor=descriptor
        source_info(self)
        # Empty dictionaries to fill later
        self.filein={}
        self.data_in={}
        self.filein_source2={}
        self.data_in_source2={}
        self.fileanncycle={}
        self.fileanncycleleap={}
        self.data_anncycle={}
        self.data_anncycleleap={}
        # 
        try:
            dummy=self.filepre
        except AttributeError:
            self.filepre=''
        #
        self.file_data_out=os.path.join(self.basedir,self.source,'std','VAR_NAME_'+str(self.level)+self.filepre+'_'+self.wildcard+'.nc')
        if self.verbose:
            print(self)        

    def __repr__(self):
        return 'CubeDiagnostics (verbose={0.verbose!r})'.format(self)

    def __str__(self):
        if self.verbose==2:
            ss=h1a+'CubeDiagnostics instance \n'+\
                'source: {0.source!s} \n'+\
                'file_data_out: {0.file_data_out!s} \n'+h1b
            return(ss.format(self))
        else:
            return 'CubeDiagnostics instance'

    def f_read_data(self,var_name,level,source=False,verbose=False):
        """Lazy read cube(s) of var_name at level for current time block.

        Add entry to the dictionary attributes self.filein and
        self.data_in.
        """
        name=var_name2long_name[var_name]
        if not source:
            source=self.source
        self.filein[var_name+'_'+str(level)]=os.path.join(self.basedir,source,'std',var_name+'_'+str(level)+self.filepre+'_'+self.wildcard+'.nc')
        self.data_in[var_name+'_'+str(level)]=iris.load(self.filein[var_name+'_'+str(level)],name)
        if verbose:
            ss=h2a+'f_read_data \n'+\
                'var_name: {0!s} \n'+\
                'filein: {1.filein!s} \n'+\
                'data_in: {1.data_in!s} \n'+h2b
            print(ss.format(var_name,self))

    def f_read_data_source2(self,var_name,level,verbose=False):
        """Lazy read cube(s) of var_name at level for current time block.

        Add entry to the dictionary attributes self.filein_source2 and
        self.data_in_source2.
        """
        name=var_name2long_name[var_name]
        self.filein_source2[var_name+'_'+str(level)]=os.path.join(self.basedir,self.source2,'std',var_name+'_'+str(level)+self.filepre+'_'+self.wildcard+'.nc')
        self.data_in_source2[var_name+'_'+str(level)]=iris.load(self.filein_source2[var_name+'_'+str(level)],name)
        if verbose:
            ss=h2a+'f_read_data_source2 \n'+\
                'var_name: {0!s} \n'+\
                'filein_source2: {1.filein_source2!s} \n'+\
                'data_in_source2: {1.data_in_source2!s} \n'+h2b
            print(ss.format(var_name,self))

    def f_read_anncycle(self,var_name,level,verbose=False):
        """Read annual cycle.

        For use generally in diagnostics where a term is decomposed into 
        annual cycle and anomaly parts, e.g.,
        PQ = Pbar Qbar + Pbar Qprime + Pprime Qbar + Pprime Qprime

        Add entry to the dictionary attributes self.fileanncycle and
        self.data_anncycle (and self.fileanncycleleap and
        self.data_anncycleleap if needed).
        """
        name=var_name2long_name[var_name]
        self.fileanncycle[var_name+'_'+str(level)]=os.path.join(self.basedir,self.source,'processed',var_name+'_'+str(level)+self.filepreanncycle+'.nc')
        self.data_anncycle[var_name+'_'+str(level)]=iris.load_cube(self.fileanncycle[var_name+'_'+str(level)],name)
        if self.calendar=='gregorian':
            self.fileanncycleleap[var_name+'_'+str(level)]=os.path.join(self.basedir,self.source,'processed',var_name+'_'+str(level)+self.filepreanncycle+'_leap.nc')
            self.data_anncycleleap[var_name+'_'+str(level)]=iris.load_cube(self.fileanncycleleap[var_name+'_'+str(level)],name)
        if verbose:
            ss=h2a+'f_read_anncycle \n'+\
                'var_name: {0!s} \n'+\
                'fileanncyle: {1.fileanncycle!s} \n'+\
                'data_anncycle: {1.data_anncycle!s} \n'
            if self.calendar=='gregorian':
                ss+='fileanncyleleap: {1.fileanncycleleap!s} \n'+\
                'data_anncycleleap: {1.data_anncycleleap!s} \n'
            ss+=h2b
            print(ss.format(var_name,self))

    def f_mld(self,method=1,zzsfc=1.0,deltatsc=1.0):
        """Calculate mixed layer depth.

        Assumes tsc (conservative temperature) has already been
        loaded, in self.data_in['tsc_all'].

        Choose <method> to calculate mixed layer depth.

        Method 1 (default).  For each profile, takes (conservative)
        temperature at the "surface" depth. (default is <zzsfc> of 1.0
        m).  Finds by linear interpolation the depth z* (zz_star) at which
        the temperature is T* (tsc_star) = <deltatsc> (default of 1.0 degC)
        below the temperature at <zzsfc>.  Saves this as attribute
        self.mltt and saves to file.

        Linear interpolation details.  For each profile, start from
        the surface and go deeper to find z_a (zz_a), the depth at
        which the temperature first falls below T*.  Temperature at
        z_a is T_a (tsc_a).  The depth immediately before z_a is then
        z_b (zz_b_, with temperature T_b (tsc_b).  T* must be in the
        range T_a < T* < T_b, and the desired z* is in the range z_b <
        z* < z_a.  From linear interpolation:

        z* = T*(z_a - z_b) - (T_b*z_a - T_a*z_b)
             ----------------------------------
                       (T_a - T_b)
                       
        
        """
        # Read in tsc for current time block and assign to tsc attribute
        self.time1,self.time2=block_times(self,verbose=self.verbose)
        time_constraint=set_time_constraint(self.time1,self.time2,calendar=self.calendar,verbose=self.verbose)
        x1=self.data_in['tsc_all'].extract(time_constraint)
        self.tsc=x1.concatenate_cube()
        # Make a copy of tsc and work on this copy
        tsc=self.tsc.copy()
        # Check that tsc has a depth axis
        dim_coord_names=[xx.var_name for xx in tsc.dim_coords]
        if 'depth' not in dim_coord_names:
            raise UserWarning('depth must be a coordinate.')
        # Transpose tsc cube such that depth is first axis
        ndim=len(dim_coord_names)
        indices=list(range(ndim))
        depth_index=dim_coord_names.index('depth')
        indices.pop(depth_index)
        indices_new=[depth_index,]+indices
        tsc.transpose(new_order=indices_new)
        # If tsc has > 2 dimensions, reshape tsc data to 2-D (nz,ngrid) numpy array
        shape1=tsc.shape
        nz=shape1[0]
        if ndim==1:
            ngrid=None
            tsc_data=tsc.data
        elif ndim==2:
            ngrid=shape1[1]
            tsc_data=tsc.data
        else:
            ngrid=shape1[1]
            for idim in range(2,ndim):
                ngrid*=shape1[idim]
            shape2=(nz,ngrid)
            tsc_data=numpy.reshape(tsc.data,shape2)
        print('ndim, shape1: {0!s}, {1!s}'.format(ndim,shape1))
        print('nz, ngrid: {0!s}, {1!s}'.format(nz,ngrid))
        #
        if method==1:
            # Method 1.
            # Extract 'surface temperature' as temperature at smallest depth.
            lev_coord=self.tsc.coord('depth')
            if lev_coord.points[0]>lev_coord.points[-1]:
                raise UserWarning('Depth coordinate must be increasing, not decreasing.')
            lev_con=iris.Constraint(depth=zzsfc)
            tsfc=tsc.extract(lev_con)
            # Create field of T* = (surface temperature minus deltatsc)
            tsc_star=tsfc.data-deltatsc
            # Find the indices of z_a, the first depth at which the temperature
            #   is less than T*
            # z_b is then the depth immediately above this.
            # T_a and T_b are the temperatures at these two depths.  They
            #   bracket T*
            indices_za=np.argmax(tsc_data<tsc_star,axis=0)
            indices_zb=indices_za-1
            tsc_a=tsc_data[indices_za,list(np.arange(ngrid))]
            tsc_b=tsc_data[indices_zb,list(np.arange(ngrid))]
            zz_a=lev_coord.points[indices_za]
            zz_b=lev_coord.points[indices_zb]
            # Linearly interpolate to find the depth z* at which the
            #   temperature is T*
            zz_star=(tsc_star*(zz_a-zz_b)-(tsc_b*zz_a-tsc_a*zz_b)) / (tsc_a-tsc_b)
            if self.verbose:
                ii=15
                print('Sample data from first grid point')
                print('deltatsc: {0!s}'.format(deltatsc))
                print('ii: {0!s}'.format(ii))
                print('tsfc.data[ii]: {0!s}'.format(tsfc.data[ii]))
                print('tsc_star[ii]: {0!s}'.format(tsc_star[ii]))
                print('indices_za[ii]: {0!s}'.format(indices_za[ii]))
                print('indices_zb[ii]: {0!s}'.format(indices_zb[ii]))
                print('tsc_a[ii]: {0!s}'.format(tsc_a[ii]))
                print('tsc_b[ii]: {0!s}'.format(tsc_b[ii]))
                print('zz_a[ii]: {0!s}'.format(zz_a[ii]))
                print('zz_b[ii]: {0!s}'.format(zz_b[ii]))
                print('zz_star[ii]: {0!s}'.format(zz_star[ii]))
        else:
            raise UserWarning('Invalid method')
        # Reshape zz_star to shape of original array (minus depth axis)
        if ndim>=2:
            zz_star.reshape(shape1[1:])
        # Create a list of two-lists, each of form [dim_coord,index]
        kdim=0
        dim_coords=[]
        for xx in tsc.dim_coords[1:]:
            dim_coords.append([xx,kdim])
            kdim+=1
        # Create iris cube of zz_star
        var_name='mltt'
        long_name=var_name2long_name[var_name]
        self.mltt=iris.cube.Cube(zz_star,units=lev_coord.units,dim_coords_and_dims=dim_coords)
        self.mltt.rename(long_name)
        self.mltt.var_name=var_name
        # Add cell method to describe calculation of mixed layer
        cm=iris.coords.CellMethod('point','depth',comments='depth where temp is temp (at depth '+str(zzsfc)+') minus '+str(deltatsc))
        self.mltt.add_cell_method(cm)
        # Save cube
        level=0
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.mltt,fileout)
        if self.archive:
            archive_file(self,fileout)

    def f_sea_water_potential_density(self,pref=0):
        """Calculate sea water potential density using Gibbs TEOS-10 routine.

        Assumes tsc (conservative temperature) and sa (absolute
        salinity) hav already been loaded, in
        self.data_in['tsc_LEVEL'] and self.data_in['sa_LEVEL'].

        Uses gsw.rho_CT(tsc,sa,pref) where pref is reference sea water
        pressure in dbar.

        Choose <method> to calculate mixed layer depth.

        Method 1 (default).  Assumes do not have explicit values of
        swp, but tsc and sa are on a grid with a depth axis (in m).
        Depth in m is sufficiently close to sea water pressure in dbar
        to use directly, when only considering shallow depths.

        Method 2.  Assumes swp has already been loaded, in
        self.data_in['swp_LEVEL'], and uses this in the calculation of
        sea water potential density.
        """
        # Read in tsc,sa for current time block and assign to tsc,sa attributes
        self.time1,self.time2=block_times(self,verbose=self.verbose)
        time_constraint=set_time_constraint(self.time1,self.time2,calendar=self.calendar,verbose=self.verbose)
        x1=self.data_in['tsc_'+str(self.level)].extract(time_constraint)
        x2=self.data_in['sa_'+str(self.level)].extract(time_constraint)
        self.tsc=x1.concatenate_cube()
        self.sa=x2.concatenate_cube()
        tsc=self.tsc.data
        sa=self.sa.data
        # Get sea water pressure
        #if method==1:
        #    # Use depth (m) as surrogate for sea water pressure (dbar)
        #    lev_coord=self.tsc.coord('depth')
        #    if lev_coord.units!='m':
        #        raise UserWarning('Depth coordinate must be in m.')
        #    swp=lev_coord.points
        #elif method==2:
        #    # Explicitly read sea water pressure data
        #    x3=self.data_in['swp_'+str(self.level)].extract(time_constraint)
        #    self.swp=x3.concatenate_cube()
        #    if self.swp.units!='dbar':
        #        raise UserWarning('Sea water pressure must be in dbar.')
        #    swp=self.swp.data
        # Calculate seawater potential density
        swpd=gsw.rho_CT_exact(sa,tsc,pref)
        # Create iris cube of swp
        var_name='swpd'
        self.swpd=create_cube(conv_float32(swpd),self.tsc,new_var_name=var_name)
        self.swpd.units='kg m-3'
        # Add cell method to describe calculation of sea water potential density
        cm=iris.coords.CellMethod('point','depth',comments='swpd calculated using gsw.rho_CT_exact in f_sea_water_potential_density: pref='+str(pref))
        self.swpd.add_cell_method(cm)
        # Save cube
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.swpd,fileout)
        if self.archive:
            archive_file(self,fileout)

    def f_specific_humidity(self):
        """Calculate specific humidity from relative humidity and temperature.

        Assumes ta (air temperature) and rhum (relative humidity) on
        pressure level self.level have already been loaded, in
        self.data_in['ta_LEVEL'] and self.data_in['rhum_LEVEL'].

        Use Clausius-Clapeyron relation (COARE algorithm) to calculate
        specific humidity.  Follows the saturation specific humidity
        computation in the COARE Fortran code v2.5b.

        Reference is
        https://woodshole.er.usgs.gov/operations/sea-mat/air_sea-html/qsat.html

        """
        # Read in ta,rhum for current time block and assign to ta,rhum attributes
        self.time1,self.time2=block_times(self,verbose=self.verbose)
        time_constraint=set_time_constraint(self.time1,self.time2,calendar=self.calendar,verbose=self.verbose)
        x1=self.data_in['ta_'+str(self.level)].extract(time_constraint)
        x2=self.data_in['rhum_'+str(self.level)].extract(time_constraint)
        self.ta=x1.concatenate_cube()
        self.rhum=x2.concatenate_cube()
        # air temperature must be in degC
        self.ta.convert_units('celsius')
        # relative humidity must be a fraction (not percentage)
        self.rhum.convert_units('1')
        # Calculate specific humidity
        ta=self.ta.data
        rhum=self.rhum.data
        pa=self.level
        # vapour pressure in hPa
        ew=6.1121*(1.0007+3.46e-6*pa)*np.exp((17.502*ta)/(240.97+ta))
        # saturated specific humidity in kg kg-1
        shum_sat=0.62197*(ew/(pa-0.378*ew))
        # specific humidity in kg kg-1
        shum=rhum*shum_sat
        # Create iris cube of shum
        var_name='shum'
        self.shum=create_cube(conv_float32(shum),self.rhum,new_var_name=var_name)
        self.shum.units='1'
        # Add cell method to describe calculation of sea water potential density
        cm=iris.coords.CellMethod('point','pressure',comments='specific humidity calculated from relative humidity and temperature')
        self.shum.add_cell_method(cm)
        # Save cube
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.shum,fileout)
        if self.archive:
            archive_file(self,fileout)
        if self.verbose:
            ta_min=self.ta.data.min()
            ta_max=self.ta.data.max()
            rhum_min=self.rhum.data.min()
            rhum_max=self.rhum.data.max()
            shum_min=self.shum.data.min()
            shum_max=self.shum.data.max()
            ss=h1a+'f_specific_humidity\n'+\
                'level: {0.level!s}\n'+\
                'ta min: {1!s}\n'+\
                'ta max: {2!s}\n'+\
                'rhum min: {3!s}\n'+\
                'rhum max: {4!s}\n'+\
                'shum min: {5!s}\n'+\
                'shum max: {6!s}\n'+h1b
            print(ss.format(self,ta_min,ta_max,rhum_min,rhum_max,shum_min,shum_max))

    def f_potential_temperature(self):
        """Calculate potential temperature from temperature and pressure.

        Assumes ta (air temperature) on pressure level self.level has
        already been loaded, in self.data_in['ta_LEVEL'].

        Uses standard definition of potential temperature (theta)

        theta = T (p_0 / p)^{R/C_p}, where

        p_0 = 10^5 Pa is standard reference pressure
        R is specific gas constant
        c_p is specific heat capacity at constant pressure.

        """
        # Read in ta for current time block and assign to ta attribute
        self.time1,self.time2=block_times(self,verbose=self.verbose)
        time_constraint=set_time_constraint(self.time1,self.time2,calendar=self.calendar,verbose=self.verbose)
        x1=self.data_in['ta_'+str(self.level)].extract(time_constraint)
        self.ta=x1.concatenate_cube()
        # air temperature must be in Kelvin
        self.ta.convert_units('Kelvin')
        # Set constants
        planet=Planet()
        rr=planet.gascon
        cp=planet.cp
        kappa=rr/cp
        p0=1e5 # reference pressure in Pascals
        pa=self.level*100 # pressure level in Pascals
        print('R,Cp,kappa: {0!s}, {1!s}, {2!s}'.format(rr,cp,kappa))
        print('p0,p: {0!s}, {1!s}'.format(p0,pa))
        factor=math.pow(p0/pa,kappa)
        print('factor: {0!s}'.format(factor))
        # Calculate potential temperature in Kelvin
        theta=self.ta.data*factor
        # Create iris cube of theta
        var_name='theta'
        self.theta=create_cube(conv_float32(theta),self.ta,new_var_name=var_name)
        self.theta.units='Kelvin'
        # Add cell method to describe calculation of sea water potential density
        cm=iris.coords.CellMethod('point','pressure',comments='potential temperature calculated from temperature and pressure')
        self.theta.add_cell_method(cm)
        # Save cube
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.theta,fileout)
        if self.archive:
            archive_file(self,fileout)

    def f_vwndptap(self):
        """Calculate v'T' from anomalous vwnd and ta.

        Assumes anomalous vwnd (northward wind) and ta (air temperature) on
        level self.level have already been loaded, in
        self.data_in['vwnd_LEVEL'] and self.data_in['ta_LEVEL'].

        vwndptap (v'T') = v' * T'
        """
        # Check that input data have had the annual cycle subtracted
        if self.filepre[:4]!='_rac':
            raise UserWarning('Input vwnd and ta must have had the annual cycle subtracted.')
        # Change the name of the output file to remove the '_rac' string
        # This is because it is understood that v'T' is calculated from anomalous data
        # The resulting v'T' maps will themselves have an annual cycle, which can be calculated
        #  using anncycle.py, along with anomalies from that annual cycle: (v'T')' !
        print('Initial file_data_out: {0.file_data_out!s}'.format(self))
        self.file_data_out=os.path.join(self.basedir,self.source,'std','VAR_NAME_'+str(self.level)+self.filepre[4:]+'_'+self.wildcard+'.nc')
        print('Final file_data_out: {0.file_data_out!s}'.format(self))
        # Read in vwnd,ta for current time block and assign to vwnd,ta attributes
        self.time1,self.time2=block_times(self,verbose=self.verbose)
        time_constraint=set_time_constraint(self.time1,self.time2,calendar=self.calendar,verbose=self.verbose)
        x1=self.data_in['vwnd_'+str(self.level)].extract(time_constraint)
        x2=self.data_in['ta_'+str(self.level)].extract(time_constraint)
        self.vwnd=x1.concatenate_cube()
        self.ta=x2.concatenate_cube()
        # Check units
        if self.vwnd.units not in ['m/s','m s-1']:
            raise UserWarning('vwnd units not valid: '+str(self.vwnd.units))
        if self.ta.units not in ['degK','degC']:
            raise UserWarning('ta units not valid: '+str(self.ta.units))
        # Calculate v'T'
        vwnd=self.vwnd.data
        ta=self.ta.data
        vwndptap=vwnd*ta
        # Create iris cube of vwndptap
        var_name='vwndptap'
        self.vwndptap=create_cube(conv_float32(vwndptap),self.vwnd,new_var_name=var_name)
        self.vwndptap.units='K m s-1'
        # Add cell method to describe calculation of sea water potential density
        cm=iris.coords.CellMethod('point','pressure',comments='northward temperature flux calculcated from anomalous northward wind and anomalous air temperature')
        self.vwndptap.add_cell_method(cm)
        # Save cube
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.vwndptap,fileout)
        if self.archive:
            archive_file(self,fileout)
        if self.verbose:
            vwnd_min=self.vwnd.data.min()
            vwnd_max=self.vwnd.data.max()
            ta_min=self.ta.data.min()
            ta_max=self.ta.data.max()
            vwndptap_min=self.vwndptap.data.min()
            vwndptap_max=self.vwndptap.data.max()
            ss=h1a+'f_vwndptap\n'+\
                'level: {0.level!s}\n'+\
                'vwnd min: {1!s}\n'+\
                'vwnd max: {2!s}\n'+\
                'ta min: {3!s}\n'+\
                'ta max: {4!s}\n'+\
                'vwndptap min: {5!s}\n'+\
                'vwndptap max: {6!s}\n'+h1b
            print(ss.format(self,vwnd_min,vwnd_max,ta_min,ta_max,vwndptap_min,vwndptap_max))

    def f_vrtbudget(self,level_below,level,level_above):
        """Calculate and save terms in vorticity budget at pressure level.

        Assumes input data has already been loaded for current time
        block, using f_read_data(), in:
        
        self.data_in['uwnd_<level_below>']
        self.data_in['uwnd_<level>']
        self.data_in['uwnd_<level_above>']
        self.data_in['vwnd_<level_below>']
        self.data_in['vwnd_<level>']
        self.data_in['vwnd_<level_above>']
        self.data_in['vrt_<level_below>']
        self.data_in['vrt_<level>']
        self.data_in['vrt_<level_above>']
        self.data_in['omega_<level>']
        self.data_in['div_<level>']

        Calculate and create attributes:

        self.south2north: True if latitude coordinate of input data
        (taken as self.data_in['uwnd_<level>']) runs from south to
        north, False otherwise.  The windspharm methods (horizontal
        gradients and calculation of ff and beta) need cubes to run
        from north to south, and convert them if they do not.  Cube
        output from windspharm running from north to south can then
        not be combined (added) with cubes running from south to
        north.  So convert cube output from windspharm to run from
        south to north if self.south2north is True.

        self.dvrtdt (d zeta/dt) by centred differences

        self.m_uwnd_dvrtdx (-u d zeta/dx)
        self.m_vwnd_dvrtdy (-v d zeta/dy)
        self.m_omega_dvrtdp (-omega d zeta/dp)

        self.m_vrt_div (-zeta D)
        self.m_ff_div (-fD)

        self.m_beta_vwnd (-beta v)

        self.m_domegadx_dvwnddp (-d omega/dx * dv/dp)
        self.domegady_duwnddp (d omega/dy * du/dp)

        self.source_dvrtdt (-u d zeta/dx -v d zeta/dy -omega d zeta/dp
                            -zeta D -fD -beta v
                           - d omega/dx * dv/dp +d omega/dy * du/dp)

        self.res_dvrtdt : d zeta/dt = source + residual
                          residual = d zeta/dt - source
                          residual is the missing source needed to create
                            the actual vorticity tendency

        """
        # Read data for current time block
        self.time1,self.time2=block_times(self,verbose=self.verbose)
        time_constraint=set_time_constraint(self.time1,self.time2,calendar=self.calendar,verbose=self.verbose)
        x1=self.data_in['uwnd_'+str(level_below)].extract(time_constraint)
        x2=self.data_in['uwnd_'+str(level)].extract(time_constraint)
        x3=self.data_in['uwnd_'+str(level_above)].extract(time_constraint)
        x4=self.data_in['vwnd_'+str(level_below)].extract(time_constraint)
        x5=self.data_in['vwnd_'+str(level)].extract(time_constraint)
        x6=self.data_in['vwnd_'+str(level_above)].extract(time_constraint)
        x7=self.data_in['vrt_'+str(level_below)].extract(time_constraint)
        x8=self.data_in['vrt_'+str(level)].extract(time_constraint)
        x9=self.data_in['vrt_'+str(level_above)].extract(time_constraint)
        x10=self.data_in['omega_'+str(level)].extract(time_constraint)
        x11=self.data_in['div_'+str(level)].extract(time_constraint)
        self.uwnd_level_below=x1.concatenate_cube()
        self.uwnd_level=x2.concatenate_cube()
        self.uwnd_level_above=x3.concatenate_cube()
        self.vwnd_level_below=x4.concatenate_cube()
        self.vwnd_level=x5.concatenate_cube()
        self.vwnd_level_above=x6.concatenate_cube()
        self.vrt_level_below=x7.concatenate_cube()
        self.vrt_level=x8.concatenate_cube()
        self.vrt_level_above=x9.concatenate_cube()
        self.omega_level=x10.concatenate_cube()
        self.div_level=x11.concatenate_cube()
        #
        # Find value of south2north
        self.south2north=f_south2north(self.uwnd_level,verbose=self.verbose)
        #
        ### Calculate dvrtdt
        lfwdfirsttime=False
        lbwdlasttime=False
        tcoord=self.vrt_level.coord('time')
        # Try to read vrt for timestep before beginning of current time block
        time_first=tcoord.units.num2date(tcoord.points[0])
        time_before=time_first-self.timedelta
        print('time_first,time_before: {0!s}, {1!s}'.format(time_first,time_before))
        time_constraint=set_time_constraint(time_before,False,calendar=self.calendar,verbose=self.verbose)
        x12=self.data_in['vrt_'+str(level)].extract(time_constraint)
        if len(x12)==0:
            print('Data at time_before not available, use forward time step for dvrtdt for first time.')
            time_constraint=set_time_constraint(time_first,False,calendar=self.calendar,verbose=self.verbose)
            x20=self.data_in['vrt_'+str(level)].extract(time_constraint)
            x21=x20.concatenate_cube()
            self.vrt_level_time_before=x21.data
            lfwdfirsttime=True
        else:
            x12a=x12.concatenate_cube()
            self.vrt_level_time_before=x12a.data
        # Try to read vrt for timestep after end of current time block
        time_last=tcoord.units.num2date(tcoord.points[-1])
        time_after=time_last+self.timedelta
        print('time_last,time_after: {0!s}, {1!s}'.format(time_last,time_after))
        time_constraint=set_time_constraint(time_after,False,calendar=self.calendar,verbose=self.verbose)
        x13=self.data_in['vrt_'+str(level)].extract(time_constraint)
        if len(x13)==0:
            print('Data at time_after not available, use backward time step for dvrtdt for last time.')
            time_constraint=set_time_constraint(time_last,False,calendar=self.calendar,verbose=self.verbose)
            x20=self.data_in['vrt_'+str(level)].extract(time_constraint)
            x21=x20.concatenate_cube()
            self.vrt_level_time_after=x21.data
            lbwdlasttime=True
        else:
            x13a=x13.concatenate_cube()
            self.vrt_level_time_after=x13a.data
        # Find index of time dimension
        dim_coord_names=[xx.var_name for xx in self.vrt_level.dim_coords]
        time_index=dim_coord_names.index('time')
        print('time_index: {0!s}'.format(time_index))
        vrt=self.vrt_level.data.copy()
        vrt_time_minus1=np.roll(vrt,1,axis=time_index)
        vrt_time_plus1=np.roll(vrt,-1,axis=time_index)
        # Don't know how to assign to an arbitrary index.
        if time_index==0:
            vrt_time_minus1[0,...]=self.vrt_level_time_before
            vrt_time_plus1[-1,...]=self.vrt_level_time_after
        else:
            raise ToDoError('Code up for time index other than 0. And also check code below.')
        # Find dvrtdt from centered differences
        deltatime_seconds=self.timedelta.total_seconds()
        print('deltatime_seconds: {0!s}'.format(deltatime_seconds))
        dvrtdt=(vrt_time_plus1-vrt_time_minus1)/(2*deltatime_seconds)
        # Use forward time step at first time if needed amd overwrite
        if lfwdfirsttime:
            dvrtdt[0,...]=dvrtdt[0,...]*2
        # Use backward time step at last time if needed and overwrite
        if lbwdlasttime:
            dvrtdt[-1,...]=dvrtdt[-1,...]*2
        # Create iris cube
        dvrtdt=conv_float32(dvrtdt)
        dvrtdt=create_cube(dvrtdt,self.vrt_level)
        var_name='dvrtdt'
        long_name=var_name2long_name[var_name]
        dvrtdt.rename(long_name) # not a standard_name
        dvrtdt.var_name=var_name
        vrt_tendency_units='s-2'
        dvrtdt.units=vrt_tendency_units
        self.dvrtdt=dvrtdt
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.dvrtdt,fileout)
        if self.archive:
            archive_file(self,fileout)
        #
        ### Calculate m_uwnd_dvrtdx and m_vwnd_dvrtdy
        ww=VectorWind(self.uwnd_level,self.vwnd_level)
        dvrtdx,dvrtdy=ww.gradient(self.vrt_level)
        if self.south2north:
            dvrtdx=lat_direction(dvrtdx,'s2n')
            dvrtdy=lat_direction(dvrtdy,'s2n')
        m_uwnd_dvrtdx=-1*self.uwnd_level*dvrtdx
        m_vwnd_dvrtdy=-1*self.vwnd_level*dvrtdy
        # m_uwnd_dvrtdx attributes
        var_name='m_uwnd_dvrtdx'
        long_name=var_name2long_name[var_name]
        m_uwnd_dvrtdx.rename(long_name) # not a standard_name
        m_uwnd_dvrtdx.var_name=var_name
        m_uwnd_dvrtdx.units=vrt_tendency_units
        self.m_uwnd_dvrtdx=m_uwnd_dvrtdx
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_uwnd_dvrtdx,fileout)
        if self.archive:
            archive_file(self,fileout)
        # m_vwnd_dvrtdy attributes
        var_name='m_vwnd_dvrtdy'
        long_name=var_name2long_name[var_name]
        m_vwnd_dvrtdy.rename(long_name) # not a standard_name
        m_vwnd_dvrtdy.var_name=var_name
        m_vwnd_dvrtdy.units=vrt_tendency_units
        self.m_vwnd_dvrtdy=m_vwnd_dvrtdy
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_vwnd_dvrtdy,fileout)
        if self.archive:
            archive_file(self,fileout)
        #
        ### Calculate m_omega_dvrtdp
        deltap=100*(level_below-level_above) # Pa
        print('deltap: {0!s}'.format(deltap))
        dvrtdp=(self.vrt_level_below-self.vrt_level_above)/deltap
        m_omega_dvrtdp=-1*self.omega_level*dvrtdp
        # Attributes
        var_name='m_omega_dvrtdp'
        long_name=var_name2long_name[var_name]
        m_omega_dvrtdp.rename(long_name) # not a standard_name
        m_omega_dvrtdp.var_name=var_name
        m_omega_dvrtdp.units=vrt_tendency_units
        self.m_omega_dvrtdp=m_omega_dvrtdp
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_omega_dvrtdp,fileout)
        if self.archive:
            archive_file(self,fileout)
        #
        ### Calculate m_vrt_div
        m_vrt_div=-1*self.vrt_level*self.div_level
        # Attributes
        var_name='m_vrt_div'
        long_name=var_name2long_name[var_name]
        m_vrt_div.rename(long_name) # not a standard_name
        m_vrt_div.var_name=var_name
        m_vrt_div.units=vrt_tendency_units
        self.m_vrt_div=m_vrt_div
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_vrt_div,fileout)
        if self.archive:
            archive_file(self,fileout)
        #
        ### Calculate m_ff_div
        ff=ww.planetaryvorticity()
        if self.south2north:
            ff=lat_direction(ff,'s2n')
        m_ff_div=-1*ff*self.div_level
        m_ff_div.data=conv_float32(m_ff_div.data)
        # Attributes
        var_name='m_ff_div'
        long_name=var_name2long_name[var_name]
        m_ff_div.rename(long_name) # not a standard_name
        m_ff_div.var_name=var_name
        m_ff_div.units=vrt_tendency_units
        self.m_ff_div=m_ff_div
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_ff_div,fileout)
        if self.archive:
            archive_file(self,fileout)
        #
        ### Calculate m_beta_vwnd
        dummy,beta=ww.gradient(ff)
        if self.south2north:
            beta=lat_direction(beta,'s2n')
        m_beta_vwnd=-1*beta*self.vwnd_level
        # Attributes
        var_name='m_beta_vwnd'
        long_name=var_name2long_name[var_name]
        m_beta_vwnd.rename(long_name) # not a standard_name
        m_beta_vwnd.var_name=var_name
        m_beta_vwnd.units=vrt_tendency_units
        self.m_beta_vwnd=m_beta_vwnd
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_beta_vwnd,fileout)
        if self.archive:
            archive_file(self,fileout)
        #
        ### Calculate m_domegadx_dvwnddp
        domegadx,domegady=ww.gradient(self.omega_level)
        if self.south2north:
            domegadx=lat_direction(domegadx,'s2n')
            domegady=lat_direction(domegady,'s2n')
        dvwnddp=(self.vwnd_level_below-self.vwnd_level_above)/deltap
        m_domegadx_dvwnddp=-1*domegadx*dvwnddp
        m_domegadx_dvwnddp.data=conv_float32(m_domegadx_dvwnddp.data)
        # Attributes
        var_name='m_domegadx_dvwnddp'
        long_name=var_name2long_name[var_name]
        m_domegadx_dvwnddp.rename(long_name) # not a standard_name
        m_domegadx_dvwnddp.var_name=var_name
        m_domegadx_dvwnddp.units=vrt_tendency_units
        self.m_domegadx_dvwnddp=m_domegadx_dvwnddp
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_domegadx_dvwnddp,fileout)
        if self.archive:
            archive_file(self,fileout)
        #
        ### Calculate domegady_duwnddp
        duwnddp=(self.uwnd_level_below-self.uwnd_level_above)/deltap
        domegady_duwnddp=domegady*duwnddp
        domegady_duwnddp.data=conv_float32(domegady_duwnddp.data)
        # Attributes
        var_name='domegady_duwnddp'
        long_name=var_name2long_name[var_name]
        domegady_duwnddp.rename(long_name) # not a standard_name
        domegady_duwnddp.var_name=var_name
        domegady_duwnddp.units=vrt_tendency_units
        self.domegady_duwnddp=domegady_duwnddp
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.domegady_duwnddp,fileout)
        if self.archive:
            archive_file(self,fileout)
        #
        ### Calculate source_dvrtdt (total source)
        source_dvrtdt=self.m_uwnd_dvrtdx+self.m_vwnd_dvrtdy+self.m_omega_dvrtdp+self.m_vrt_div+self.m_ff_div+self.m_beta_vwnd+self.m_domegadx_dvwnddp+self.domegady_duwnddp
        source_dvrtdt.data=conv_float32(source_dvrtdt.data)
        # Attributes
        var_name='source_dvrtdt'
        long_name=var_name2long_name[var_name]
        source_dvrtdt.rename(long_name) # not a standard_name
        source_dvrtdt.var_name=var_name
        source_dvrtdt.units=vrt_tendency_units
        self.source_dvrtdt=source_dvrtdt
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.source_dvrtdt,fileout)
        if self.archive:
            archive_file(self,fileout)
        #
        ### Calculate res_dvrtdt (residual)
        res_dvrtdt=self.dvrtdt-self.source_dvrtdt
        res_dvrtdt.data=conv_float32(res_dvrtdt.data)
        # Attributes
        var_name='res_dvrtdt'
        long_name=var_name2long_name[var_name]
        res_dvrtdt.rename(long_name) # not a standard_name
        res_dvrtdt.var_name=var_name
        res_dvrtdt.units=vrt_tendency_units
        self.res_dvrtdt=res_dvrtdt
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.res_dvrtdt,fileout)
        if self.archive:
            archive_file(self,fileout)

    def f_m_vrt_div(self,level):
        """Calculate and save m_vrt_div term only in vorticity budget at pressure level.

        This is a subset of the code in f_vrtbudget, just to calculate the single term 
        m_vrt_div. For full documentation, see f_vrtbudget.

        """
        # Read data for current time block
        self.time1,self.time2=block_times(self,verbose=self.verbose)
        time_constraint=set_time_constraint(self.time1,self.time2,calendar=self.calendar,verbose=self.verbose)
        x8=self.data_in['vrt_'+str(level)].extract(time_constraint)
        x11=self.data_in['div_'+str(level)].extract(time_constraint)
        self.vrt_level=x8.concatenate_cube()
        self.div_level=x11.concatenate_cube()
        #
        ### Calculate m_vrt_div
        m_vrt_div=-1*self.vrt_level*self.div_level
        # Attributes
        var_name='m_vrt_div'
        long_name=var_name2long_name[var_name]
        m_vrt_div.rename(long_name) # not a standard_name
        m_vrt_div.var_name=var_name
        vrt_tendency_units='s-2'
        m_vrt_div.units=vrt_tendency_units
        self.m_vrt_div=m_vrt_div
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_vrt_div,fileout)
        if self.archive:
            archive_file(self,fileout)

    def f_m_vrt_div_annpert(self,level):
        """Calculate and save m_vrt_div terms decomposed into anncycle and perturbation parts.

        Both vorticity (vrt) and divergence (div) terms have been
        previously separated into annual cycle and perturbation parts:

        vrt = vrtbar + vrtprime
        div = divbar + divprime

        where bar refers to annual cycle part and prime refers to anomaly from annual cycle.

        m_vrt_div = -1*vrt*div
                  = -1*(vrtbar+vrtprime)*(divbar+divprime)
                  = -1*vrtbar*divbar -1*vrtbar*divprime -1*vrtprime*divbar -1*vrtprime*vrtprime
                  =  m_vrtbar_divbar +m_vrtbar_divprime =m_vrtprime_divbar +m_vrtprime_vrtprime


        Calculate attributes:

        m_vrtbar_divbar
        m_vrtbar_divprime
        m_vrtprime_divbar
        m_vrtprime_divprime
        """
        # Read vrtprime and divprime data for current time block
        if self.filepre!='_rac':
            raise UserWarning('Must read data with annual cycle removed, for perturbation.')
        self.time1,self.time2=block_times(self,verbose=self.verbose)
        time_constraint=set_time_constraint(self.time1,self.time2,calendar=self.calendar,verbose=self.verbose)
        x1=self.data_in['vrt_'+str(level)].extract(time_constraint)
        x2=self.data_in['div_'+str(level)].extract(time_constraint)
        vrtprime=x1.concatenate_cube()
        divprime=x2.concatenate_cube()
        #
        # Set annual cycles for this year (regular or leap year) and extract data for current time block
        vrtbar=self.data_anncycle['vrt_'+str(level)]
        divbar=self.data_anncycle['div_'+str(level)]
        if divmod(self.year,4)[1]==0 and self.calendar=='gregorian':
            vrtbar=self.data_anncycleleap['vrt_'+str(level)]
            divbar=self.data_anncycleleap['div_'+str(level)]
        tcoord_anncycle=vrtbar.coord('time')
        anncycle_year=tcoord_anncycle.units.num2date(tcoord_anncycle.points[0]).year
        xx=self.time1
        x1=datetime.datetime(anncycle_year,xx.month,xx.day,xx.hour,xx.minute)
        xx=self.time2
        x2=datetime.datetime(anncycle_year,xx.month,xx.day,xx.hour,xx.minute)
        print('anncycle_year,x1,x2: {0!s}, {1!s}, {2!s}'.format(anncycle_year,x1,x2))
        time_anncycle_constraint=set_time_constraint(x1,x2,calendar=self.calendar)
        vrtbar=vrtbar.extract(time_anncycle_constraint)
        divbar=divbar.extract(time_anncycle_constraint)
        #
        ### Calculate m_vrtbar_divbar
        m_vrtbar_divbar=-1*vrtbar*divbar
        # Reset time axis to current year
        m_vrtbar_divbar=create_cube(m_vrtbar_divbar.data,vrtprime)
        # Attributes
        var_name='m_vrtbar_divbar'
        long_name=var_name2long_name[var_name]
        m_vrtbar_divbar.rename(long_name) # not a standard_name
        m_vrtbar_divbar.var_name=var_name
        vrt_tendency_units='s-2'
        m_vrtbar_divbar.units=vrt_tendency_units
        self.m_vrtbar_divbar=m_vrtbar_divbar
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        fileout=fileout.replace(self.filepre,'')
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_vrtbar_divbar,fileout)
        if self.archive:
            archive_file(self,fileout)
        #
        ### Calculate m_vrtbar_divprime
        m_vrtbar_divprime=-1*vrtbar.data*divprime.data
        # Set time axis to current year
        m_vrtbar_divprime=create_cube(m_vrtbar_divprime,vrtprime)
        # Attributes
        var_name='m_vrtbar_divprime'
        long_name=var_name2long_name[var_name]
        m_vrtbar_divprime.rename(long_name) # not a standard_name
        m_vrtbar_divprime.var_name=var_name
        m_vrtbar_divprime.units=vrt_tendency_units
        self.m_vrtbar_divprime=m_vrtbar_divprime
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        fileout=fileout.replace(self.filepre,'')
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_vrtbar_divprime,fileout)
        if self.archive:
            archive_file(self,fileout)
        #
        ### Calculate m_vrtprime_divbar
        m_vrtprime_divbar=-1*vrtprime.data*divbar.data
        # Set time axis to current year
        m_vrtprime_divbar=create_cube(m_vrtprime_divbar,vrtprime)
        # Attributes
        var_name='m_vrtprime_divbar'
        long_name=var_name2long_name[var_name]
        m_vrtprime_divbar.rename(long_name) # not a standard_name
        m_vrtprime_divbar.var_name=var_name
        m_vrtprime_divbar.units=vrt_tendency_units
        self.m_vrtprime_divbar=m_vrtprime_divbar
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        fileout=fileout.replace(self.filepre,'')
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_vrtprime_divbar,fileout)
        if self.archive:
            archive_file(self,fileout)
        #
        ### Calculate m_vrtprime_divprime
        m_vrtprime_divprime=-1*vrtprime*divprime
        # Attributes
        var_name='m_vrtprime_divprime'
        long_name=var_name2long_name[var_name]
        m_vrtprime_divprime.rename(long_name) # not a standard_name
        m_vrtprime_divprime.var_name=var_name
        m_vrtprime_divprime.units=vrt_tendency_units
        self.m_vrtprime_divprime=m_vrtprime_divprime
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        fileout=fileout.replace(self.filepre,'')
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_vrtprime_divprime,fileout)
        if self.archive:
            archive_file(self,fileout)

    def f_m_vrt_div_pertremainder(self,level):
        """Calculate and save m_vrt_div terms decomposed into remainder and perturbation parts.

        This is basically the same as f_m_vrt_div_annpert, except that
        the prime terms are now from some kind of perturbation or
        filtered (e.g., CCKW filtered, EK1) input vrt and div, and the
        bar terms are from the remainder (e.g., total vrt minus the
        filtered vrt, calculated in subtract.py, e.g., NEK1).

        Both vorticity (vrt) and divergence (div) terms have been
        previously separated into remainder and perturbation parts
        (remainder is total minus perturbation, calculated using subtract.py):

        vrt = vrtbar + vrtprime
        div = divbar + divprime

        where bar refers to remainder part and prime refers to perturbation part.

        m_vrt_div = -1*vrt*div
                  = -1*(vrtbar+vrtprime)*(divbar+divprime)
                  = -1*vrtbar*divbar -1*vrtbar*divprime -1*vrtprime*divbar -1*vrtprime*vrtprime
                  =  m_vrtbar_divbar +m_vrtbar_divprime =m_vrtprime_divbar +m_vrtprime_vrtprime


        Calculate attributes:

        m_vrtbar_divbar
        m_vrtbar_divprime
        m_vrtprime_divbar
        m_vrtprime_divprime

        """
        # Read vrtprime and divprime data for current time block
        self.time1,self.time2=block_times(self,verbose=self.verbose)
        time_constraint=set_time_constraint(self.time1,self.time2,calendar=self.calendar,verbose=self.verbose)
        x1=self.data_in['vrt_'+str(level)].extract(time_constraint)
        x2=self.data_in['div_'+str(level)].extract(time_constraint)
        vrtprime=x1.concatenate_cube()
        divprime=x2.concatenate_cube()
        #
        # Read vrtbar and divbar data for current time block
        x1=self.data_in_source2['vrt_'+str(level)].extract(time_constraint)
        x2=self.data_in_source2['div_'+str(level)].extract(time_constraint)
        vrtbar=x1.concatenate_cube()
        divbar=x2.concatenate_cube()
        #
        # Following code is identical to that in f_m_vrt_div_annpert
        #
        ### Calculate m_vrtbar_divbar
        m_vrtbar_divbar=-1*vrtbar*divbar
        m_vrtbar_divbar=create_cube(m_vrtbar_divbar.data,vrtprime)
        # Attributes
        var_name='m_vrtbar_divbar'
        long_name=var_name2long_name[var_name]
        m_vrtbar_divbar.rename(long_name) # not a standard_name
        m_vrtbar_divbar.var_name=var_name
        vrt_tendency_units='s-2'
        m_vrtbar_divbar.units=vrt_tendency_units
        self.m_vrtbar_divbar=m_vrtbar_divbar
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        fileout=fileout.replace(self.filepre,'')
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_vrtbar_divbar,fileout)
        if self.archive:
            archive_file(self,fileout)
        #
        ### Calculate m_vrtbar_divprime
        m_vrtbar_divprime=-1*vrtbar.data*divprime.data
        m_vrtbar_divprime=create_cube(m_vrtbar_divprime,vrtprime)
        # Attributes
        var_name='m_vrtbar_divprime'
        long_name=var_name2long_name[var_name]
        m_vrtbar_divprime.rename(long_name) # not a standard_name
        m_vrtbar_divprime.var_name=var_name
        m_vrtbar_divprime.units=vrt_tendency_units
        self.m_vrtbar_divprime=m_vrtbar_divprime
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        fileout=fileout.replace(self.filepre,'')
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_vrtbar_divprime,fileout)
        if self.archive:
            archive_file(self,fileout)
        #
        ### Calculate m_vrtprime_divbar
        m_vrtprime_divbar=-1*vrtprime.data*divbar.data
        m_vrtprime_divbar=create_cube(m_vrtprime_divbar,vrtprime)
        # Attributes
        var_name='m_vrtprime_divbar'
        long_name=var_name2long_name[var_name]
        m_vrtprime_divbar.rename(long_name) # not a standard_name
        m_vrtprime_divbar.var_name=var_name
        m_vrtprime_divbar.units=vrt_tendency_units
        self.m_vrtprime_divbar=m_vrtprime_divbar
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        fileout=fileout.replace(self.filepre,'')
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_vrtprime_divbar,fileout)
        if self.archive:
            archive_file(self,fileout)
        #
        ### Calculate m_vrtprime_divprime
        m_vrtprime_divprime=-1*vrtprime*divprime
        # Attributes
        var_name='m_vrtprime_divprime'
        long_name=var_name2long_name[var_name]
        m_vrtprime_divprime.rename(long_name) # not a standard_name
        m_vrtprime_divprime.var_name=var_name
        m_vrtprime_divprime.units=vrt_tendency_units
        self.m_vrtprime_divprime=m_vrtprime_divprime
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        fileout=fileout.replace(self.filepre,'')
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_vrtprime_divprime,fileout)
        if self.archive:
            archive_file(self,fileout)

    def f_m_uwnd_dvrtdx_annpert(self,level):
        """Calculate and save m_uwnd_dvrtdx terms decomposed into anncycle and perturbation parts.

        Both uwnd and vorticity (vrt) terms have been
        previously separated into annual cycle and perturbation parts:

        Behaves similarly to f_m_vrt_div_annpert()

        Calculate attributes:

        m_uwndbar_dvrtdxbar
        m_uwndbar_dvrtdxprime
        m_uwndprime_dvrtdxbar
        m_uwndprime_dvrtdxprime
        """
        # Read uwndprime and vrtprime data for current time block
        if self.filepre!='_rac':
            raise UserWarning('Must read data with annual cycle removed, for perturbation.')
        self.time1,self.time2=block_times(self,verbose=self.verbose)
        time_constraint=set_time_constraint(self.time1,self.time2,calendar=self.calendar,verbose=self.verbose)
        x1=self.data_in['uwnd_'+str(level)].extract(time_constraint)
        x2=self.data_in['vrt_'+str(level)].extract(time_constraint)
        uwndprime=x1.concatenate_cube()
        vrtprime=x2.concatenate_cube()
        #
        # Set annual cycles for this year (regular or leap year) and extract data for current time block
        uwndbar=self.data_anncycle['uwnd_'+str(level)]
        vrtbar=self.data_anncycle['vrt_'+str(level)]
        if divmod(self.year,4)[1]==0 and self.calendar=='gregorian':
            uwndbar=self.data_anncycleleap['uwnd_'+str(level)]
            vrtbar=self.data_anncycleleap['vrt_'+str(level)]
        tcoord_anncycle=uwndbar.coord('time')
        anncycle_year=tcoord_anncycle.units.num2date(tcoord_anncycle.points[0]).year
        xx=self.time1
        x1=datetime.datetime(anncycle_year,xx.month,xx.day,xx.hour,xx.minute)
        xx=self.time2
        x2=datetime.datetime(anncycle_year,xx.month,xx.day,xx.hour,xx.minute)
        print('anncycle_year,x1,x2: {0!s}, {1!s}, {2!s}'.format(anncycle_year,x1,x2))
        time_anncycle_constraint=set_time_constraint(x1,x2,calendar=self.calendar)
        uwndbar=uwndbar.extract(time_anncycle_constraint)
        vrtbar=vrtbar.extract(time_anncycle_constraint)
        #
        # Find value of south2north
        self.south2north=f_south2north(uwndbar,verbose=self.verbose)
        #
        # Calculate dvrtdxbar
        # ww is dummy VectorWind instance using available data: uwndbar and vrtbar!!
        ww=VectorWind(uwndbar,vrtbar)
        dvrtdxbar,dvrtdybar=ww.gradient(vrtbar)
        if self.south2north:
            dvrtdxbar=lat_direction(dvrtdxbar,'s2n')
            dvrtdybar=lat_direction(dvrtdybar,'s2n')
        #
        # Calculate dvrtdxprime
        # ww is dummy VectorWind instance using available data: uwndprime and vrtprime!!
        ww=VectorWind(uwndprime,vrtprime)
        dvrtdxprime,dvrtdyprime=ww.gradient(vrtprime)
        if self.south2north:
            dvrtdxprime=lat_direction(dvrtdxprime,'s2n')
            dvrtdyprime=lat_direction(dvrtdyprime,'s2n')
        #
        ### Calculate m_uwndbar_dvrtdxbar
        m_uwndbar_dvrtdxbar=-1*uwndbar*dvrtdxbar
        # Reset time axis to current year
        m_uwndbar_dvrtdxbar=create_cube(m_uwndbar_dvrtdxbar.data,vrtprime)
        # Attributes
        var_name='m_uwndbar_dvrtdxbar'
        long_name=var_name2long_name[var_name]
        m_uwndbar_dvrtdxbar.rename(long_name) # not a standard_name
        m_uwndbar_dvrtdxbar.var_name=var_name
        vrt_tendency_units='s-2'
        m_uwndbar_dvrtdxbar.units=vrt_tendency_units
        self.m_uwndbar_dvrtdxbar=m_uwndbar_dvrtdxbar
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        fileout=fileout.replace(self.filepre,'')
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_uwndbar_dvrtdxbar,fileout)
        if self.archive:
            archive_file(self,fileout)
        #
        ### Calculate m_uwndbar_dvrtdxprime
        m_uwndbar_dvrtdxprime=-1*uwndbar.data*dvrtdxprime.data
        # Set time axis to current year
        m_uwndbar_dvrtdxprime=create_cube(m_uwndbar_dvrtdxprime,vrtprime)
        # Attributes
        var_name='m_uwndbar_dvrtdxprime'
        long_name=var_name2long_name[var_name]
        m_uwndbar_dvrtdxprime.rename(long_name) # not a standard_name
        m_uwndbar_dvrtdxprime.var_name=var_name
        m_uwndbar_dvrtdxprime.units=vrt_tendency_units
        self.m_uwndbar_dvrtdxprime=m_uwndbar_dvrtdxprime
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        fileout=fileout.replace(self.filepre,'')
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_uwndbar_dvrtdxprime,fileout)
        if self.archive:
            archive_file(self,fileout)
        #
        ### Calculate m_uwndprime_dvrtdxbar
        m_uwndprime_dvrtdxbar=-1*uwndprime.data*dvrtdxbar.data
        # Set time axis to current year
        m_uwndprime_dvrtdxbar=create_cube(m_uwndprime_dvrtdxbar,vrtprime)
        # Attributes
        var_name='m_uwndprime_dvrtdxbar'
        long_name=var_name2long_name[var_name]
        m_uwndprime_dvrtdxbar.rename(long_name) # not a standard_name
        m_uwndprime_dvrtdxbar.var_name=var_name
        m_uwndprime_dvrtdxbar.units=vrt_tendency_units
        self.m_uwndprime_dvrtdxbar=m_uwndprime_dvrtdxbar
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        fileout=fileout.replace(self.filepre,'')
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_uwndprime_dvrtdxbar,fileout)
        if self.archive:
            archive_file(self,fileout)
        #
        ### Calculate m_uwndprime_dvrtdxprime
        m_uwndprime_dvrtdxprime=-1*uwndprime*dvrtdxprime
        # Attributes
        var_name='m_uwndprime_dvrtdxprime'
        long_name=var_name2long_name[var_name]
        m_uwndprime_dvrtdxprime.rename(long_name) # not a standard_name
        m_uwndprime_dvrtdxprime.var_name=var_name
        m_uwndprime_dvrtdxprime.units=vrt_tendency_units
        self.m_uwndprime_dvrtdxprime=m_uwndprime_dvrtdxprime
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        fileout=fileout.replace(self.filepre,'')
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_uwndprime_dvrtdxprime,fileout)
        if self.archive:
            archive_file(self,fileout)

    def f_m_uwnd_dvrtdx_pertremainder(self,level):
        """Calculate and save m_uwnd_dvrtdx terms decomposed into remainder and perturbation parts.

        Behaves similarly to f_m_vrt_div_pertremainder()

        Calculate attributes:

        m_uwndbar_dvrtdxbar
        m_uwndbar_dvrtdxprime
        m_uwndprime_dvrtdxbar
        m_uwndprime_dvrtdxprime
        """
        # Read uwndprime and vrtprime data for current time block
        self.time1,self.time2=block_times(self,verbose=self.verbose)
        time_constraint=set_time_constraint(self.time1,self.time2,calendar=self.calendar,verbose=self.verbose)
        x1=self.data_in['uwnd_'+str(level)].extract(time_constraint)
        x2=self.data_in['vrt_'+str(level)].extract(time_constraint)
        uwndprime=x1.concatenate_cube()
        vrtprime=x2.concatenate_cube()
        #
        # Read uwndbar and vrtbar data for current time block
        x1=self.data_in_source2['uwnd_'+str(level)].extract(time_constraint)
        x2=self.data_in_source2['vrt_'+str(level)].extract(time_constraint)
        uwndbar=x1.concatenate_cube()
        vrtbar=x2.concatenate_cube()
        #
        # Following code is identical to that in f_m_uwnd_dvrtdx_annpert
        #
        # Find value of south2north
        self.south2north=f_south2north(uwndbar,verbose=self.verbose)
        #
        # Calculate dvrtdxbar
        # ww is dummy VectorWind instance using available data: uwndbar and vrtbar!!
        ww=VectorWind(uwndbar,vrtbar)
        dvrtdxbar,dvrtdybar=ww.gradient(vrtbar)
        if self.south2north:
            dvrtdxbar=lat_direction(dvrtdxbar,'s2n')
            dvrtdybar=lat_direction(dvrtdybar,'s2n')
        #
        # Calculate dvrtdxprime
        # ww is dummy VectorWind instance using available data: uwndprime and vrtprime!!
        ww=VectorWind(uwndprime,vrtprime)
        dvrtdxprime,dvrtdyprime=ww.gradient(vrtprime)
        if self.south2north:
            dvrtdxprime=lat_direction(dvrtdxprime,'s2n')
            dvrtdyprime=lat_direction(dvrtdyprime,'s2n')
        #
        ### Calculate m_uwndbar_dvrtdxbar
        m_uwndbar_dvrtdxbar=-1*uwndbar*dvrtdxbar
        m_uwndbar_dvrtdxbar=create_cube(m_uwndbar_dvrtdxbar.data,vrtprime)
        # Attributes
        var_name='m_uwndbar_dvrtdxbar'
        long_name=var_name2long_name[var_name]
        m_uwndbar_dvrtdxbar.rename(long_name) # not a standard_name
        m_uwndbar_dvrtdxbar.var_name=var_name
        vrt_tendency_units='s-2'
        m_uwndbar_dvrtdxbar.units=vrt_tendency_units
        self.m_uwndbar_dvrtdxbar=m_uwndbar_dvrtdxbar
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        fileout=fileout.replace(self.filepre,'')
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_uwndbar_dvrtdxbar,fileout)
        if self.archive:
            archive_file(self,fileout)
        #
        ### Calculate m_uwndbar_dvrtdxprime
        m_uwndbar_dvrtdxprime=-1*uwndbar.data*dvrtdxprime.data
        m_uwndbar_dvrtdxprime=create_cube(m_uwndbar_dvrtdxprime,vrtprime)
        # Attributes
        var_name='m_uwndbar_dvrtdxprime'
        long_name=var_name2long_name[var_name]
        m_uwndbar_dvrtdxprime.rename(long_name) # not a standard_name
        m_uwndbar_dvrtdxprime.var_name=var_name
        m_uwndbar_dvrtdxprime.units=vrt_tendency_units
        self.m_uwndbar_dvrtdxprime=m_uwndbar_dvrtdxprime
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        fileout=fileout.replace(self.filepre,'')
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_uwndbar_dvrtdxprime,fileout)
        if self.archive:
            archive_file(self,fileout)
        #
        ### Calculate m_uwndprime_dvrtdxbar
        m_uwndprime_dvrtdxbar=-1*uwndprime.data*dvrtdxbar.data
        m_uwndprime_dvrtdxbar=create_cube(m_uwndprime_dvrtdxbar,vrtprime)
        # Attributes
        var_name='m_uwndprime_dvrtdxbar'
        long_name=var_name2long_name[var_name]
        m_uwndprime_dvrtdxbar.rename(long_name) # not a standard_name
        m_uwndprime_dvrtdxbar.var_name=var_name
        m_uwndprime_dvrtdxbar.units=vrt_tendency_units
        self.m_uwndprime_dvrtdxbar=m_uwndprime_dvrtdxbar
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        fileout=fileout.replace(self.filepre,'')
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_uwndprime_dvrtdxbar,fileout)
        if self.archive:
            archive_file(self,fileout)
        #
        ### Calculate m_uwndprime_dvrtdxprime
        m_uwndprime_dvrtdxprime=-1*uwndprime*dvrtdxprime
        # Attributes
        var_name='m_uwndprime_dvrtdxprime'
        long_name=var_name2long_name[var_name]
        m_uwndprime_dvrtdxprime.rename(long_name) # not a standard_name
        m_uwndprime_dvrtdxprime.var_name=var_name
        m_uwndprime_dvrtdxprime.units=vrt_tendency_units
        self.m_uwndprime_dvrtdxprime=m_uwndprime_dvrtdxprime
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        fileout=fileout.replace(self.filepre,'')
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_uwndprime_dvrtdxprime,fileout)
        if self.archive:
            archive_file(self,fileout)

    def f_m_vwnd_dvrtdy_annpert(self,level):
        """Calculate and save m_vwnd_dvrtdy terms decomposed into anncycle and perturbation parts.

        Both vwnd and vorticity (vrt) terms have been
        previously separated into annual cycle and perturbation parts:

        Behaves similarly to f_m_uwnd_dvrtdx_annpert()

        Calculate attributes:

        m_vwndbar_dvrtdybar
        m_vwndbar_dvrtdyprime
        m_vwndprime_dvrtdybar
        m_vwndprime_dvrtdyprime
        """
        # Read vwndprime and vrtprime data for current time block
        if self.filepre!='_rac':
            raise UserWarning('Must read data with annual cycle removed, for perturbation.')
        self.time1,self.time2=block_times(self,verbose=self.verbose)
        time_constraint=set_time_constraint(self.time1,self.time2,calendar=self.calendar,verbose=self.verbose)
        x1=self.data_in['vwnd_'+str(level)].extract(time_constraint)
        x2=self.data_in['vrt_'+str(level)].extract(time_constraint)
        vwndprime=x1.concatenate_cube()
        vrtprime=x2.concatenate_cube()
        #
        # Set annual cycles for this year (regular or leap year) and extract data for current time block
        vwndbar=self.data_anncycle['vwnd_'+str(level)]
        vrtbar=self.data_anncycle['vrt_'+str(level)]
        if divmod(self.year,4)[1]==0 and self.calendar=='gregorian':
            vwndbar=self.data_anncycleleap['vwnd_'+str(level)]
            vrtbar=self.data_anncycleleap['vrt_'+str(level)]
        tcoord_anncycle=vwndbar.coord('time')
        anncycle_year=tcoord_anncycle.units.num2date(tcoord_anncycle.points[0]).year
        xx=self.time1
        x1=datetime.datetime(anncycle_year,xx.month,xx.day,xx.hour,xx.minute)
        xx=self.time2
        x2=datetime.datetime(anncycle_year,xx.month,xx.day,xx.hour,xx.minute)
        print('anncycle_year,x1,x2: {0!s}, {1!s}, {2!s}'.format(anncycle_year,x1,x2))
        time_anncycle_constraint=set_time_constraint(x1,x2,calendar=self.calendar)
        vwndbar=vwndbar.extract(time_anncycle_constraint)
        vrtbar=vrtbar.extract(time_anncycle_constraint)
        #
        # Find value of south2north
        self.south2north=f_south2north(vwndbar,verbose=self.verbose)
        #
        # Calculate dvrtdybar
        # ww is dummy VectorWind instance using available data: vwndbar and vrtbar!!
        ww=VectorWind(vwndbar,vrtbar)
        dvrtdxbar,dvrtdybar=ww.gradient(vrtbar)
        if self.south2north:
            dvrtdxbar=lat_direction(dvrtdxbar,'s2n')
            dvrtdybar=lat_direction(dvrtdybar,'s2n')
        #
        # Calculate dvrtdyprime
        # ww is dummy VectorWind instance using available data: vwndprime and vrtprime!!
        ww=VectorWind(vwndprime,vrtprime)
        dvrtdxprime,dvrtdyprime=ww.gradient(vrtprime)
        if self.south2north:
            dvrtdxprime=lat_direction(dvrtdxprime,'s2n')
            dvrtdyprime=lat_direction(dvrtdyprime,'s2n')
        #
        ### Calculate m_vwndbar_dvrtdybar
        m_vwndbar_dvrtdybar=-1*vwndbar*dvrtdybar
        # Reset time axis to current year
        m_vwndbar_dvrtdybar=create_cube(m_vwndbar_dvrtdybar.data,vrtprime)
        # Attributes
        var_name='m_vwndbar_dvrtdybar'
        long_name=var_name2long_name[var_name]
        m_vwndbar_dvrtdybar.rename(long_name) # not a standard_name
        m_vwndbar_dvrtdybar.var_name=var_name
        vrt_tendency_units='s-2'
        m_vwndbar_dvrtdybar.units=vrt_tendency_units
        self.m_vwndbar_dvrtdybar=m_vwndbar_dvrtdybar
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        fileout=fileout.replace(self.filepre,'')
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_vwndbar_dvrtdybar,fileout)
        if self.archive:
            archive_file(self,fileout)
        #
        ### Calculate m_vwndbar_dvrtdyprime
        m_vwndbar_dvrtdyprime=-1*vwndbar.data*dvrtdyprime.data
        # Set time axis to current year
        m_vwndbar_dvrtdyprime=create_cube(m_vwndbar_dvrtdyprime,vrtprime)
        # Attributes
        var_name='m_vwndbar_dvrtdyprime'
        long_name=var_name2long_name[var_name]
        m_vwndbar_dvrtdyprime.rename(long_name) # not a standard_name
        m_vwndbar_dvrtdyprime.var_name=var_name
        m_vwndbar_dvrtdyprime.units=vrt_tendency_units
        self.m_vwndbar_dvrtdyprime=m_vwndbar_dvrtdyprime
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        fileout=fileout.replace(self.filepre,'')
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_vwndbar_dvrtdyprime,fileout)
        if self.archive:
            archive_file(self,fileout)
        #
        ### Calculate m_vwndprime_dvrtdybar
        m_vwndprime_dvrtdybar=-1*vwndprime.data*dvrtdybar.data
        # Set time axis to current year
        m_vwndprime_dvrtdybar=create_cube(m_vwndprime_dvrtdybar,vrtprime)
        # Attributes
        var_name='m_vwndprime_dvrtdybar'
        long_name=var_name2long_name[var_name]
        m_vwndprime_dvrtdybar.rename(long_name) # not a standard_name
        m_vwndprime_dvrtdybar.var_name=var_name
        m_vwndprime_dvrtdybar.units=vrt_tendency_units
        self.m_vwndprime_dvrtdybar=m_vwndprime_dvrtdybar
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        fileout=fileout.replace(self.filepre,'')
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_vwndprime_dvrtdybar,fileout)
        if self.archive:
            archive_file(self,fileout)
        #
        ### Calculate m_vwndprime_dvrtdyprime
        m_vwndprime_dvrtdyprime=-1*vwndprime*dvrtdyprime
        # Attributes
        var_name='m_vwndprime_dvrtdyprime'
        long_name=var_name2long_name[var_name]
        m_vwndprime_dvrtdyprime.rename(long_name) # not a standard_name
        m_vwndprime_dvrtdyprime.var_name=var_name
        m_vwndprime_dvrtdyprime.units=vrt_tendency_units
        self.m_vwndprime_dvrtdyprime=m_vwndprime_dvrtdyprime
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        fileout=fileout.replace(self.filepre,'')
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_vwndprime_dvrtdyprime,fileout)
        if self.archive:
            archive_file(self,fileout)

    def f_m_vwnd_dvrtdy_pertremainder(self,level):
        """Calculate and save m_vwnd_dvrtdy terms decomposed into remainder and perturbation parts.

        Behaves similarly to f_m_uwnd_dvrtdx_pertremainder()

        Calculate attributes:

        m_vwndbar_dvrtdybar
        m_vwndbar_dvrtdyprime
        m_vwndprime_dvrtdybar
        m_vwndprime_dvrtdyprime
        """
        # Read vwndprime and vrtprime data for current time block
        self.time1,self.time2=block_times(self,verbose=self.verbose)
        time_constraint=set_time_constraint(self.time1,self.time2,calendar=self.calendar,verbose=self.verbose)
        x1=self.data_in['vwnd_'+str(level)].extract(time_constraint)
        x2=self.data_in['vrt_'+str(level)].extract(time_constraint)
        vwndprime=x1.concatenate_cube()
        vrtprime=x2.concatenate_cube()
        #
        # Read vwndbar and vrtbar data for current time block
        x1=self.data_in_source2['vwnd_'+str(level)].extract(time_constraint)
        x2=self.data_in_source2['vrt_'+str(level)].extract(time_constraint)
        vwndbar=x1.concatenate_cube()
        vrtbar=x2.concatenate_cube()
        #
        # Following code is identical to that in f_m_vwnd_dvrtdy_annpert
        #
        # Find value of south2north
        self.south2north=f_south2north(vwndbar,verbose=self.verbose)
        #
        # Calculate dvrtdybar
        # ww is dummy VectorWind instance using available data: vwndbar and vrtbar!!
        ww=VectorWind(vwndbar,vrtbar)
        dvrtdxbar,dvrtdybar=ww.gradient(vrtbar)
        if self.south2north:
            dvrtdxbar=lat_direction(dvrtdxbar,'s2n')
            dvrtdybar=lat_direction(dvrtdybar,'s2n')
        #
        # Calculate dvrtdyprime
        # ww is dummy VectorWind instance using available data: vwndprime and vrtprime!!
        ww=VectorWind(vwndprime,vrtprime)
        dvrtdxprime,dvrtdyprime=ww.gradient(vrtprime)
        if self.south2north:
            dvrtdxprime=lat_direction(dvrtdxprime,'s2n')
            dvrtdyprime=lat_direction(dvrtdyprime,'s2n')
        #
        ### Calculate m_vwndbar_dvrtdybar
        m_vwndbar_dvrtdybar=-1*vwndbar*dvrtdybar
        m_vwndbar_dvrtdybar=create_cube(m_vwndbar_dvrtdybar.data,vrtprime)
        # Attributes
        var_name='m_vwndbar_dvrtdybar'
        long_name=var_name2long_name[var_name]
        m_vwndbar_dvrtdybar.rename(long_name) # not a standard_name
        m_vwndbar_dvrtdybar.var_name=var_name
        vrt_tendency_units='s-2'
        m_vwndbar_dvrtdybar.units=vrt_tendency_units
        self.m_vwndbar_dvrtdybar=m_vwndbar_dvrtdybar
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        fileout=fileout.replace(self.filepre,'')
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_vwndbar_dvrtdybar,fileout)
        if self.archive:
            archive_file(self,fileout)
        #
        ### Calculate m_vwndbar_dvrtdyprime
        m_vwndbar_dvrtdyprime=-1*vwndbar.data*dvrtdyprime.data
        m_vwndbar_dvrtdyprime=create_cube(m_vwndbar_dvrtdyprime,vrtprime)
        # Attributes
        var_name='m_vwndbar_dvrtdyprime'
        long_name=var_name2long_name[var_name]
        m_vwndbar_dvrtdyprime.rename(long_name) # not a standard_name
        m_vwndbar_dvrtdyprime.var_name=var_name
        m_vwndbar_dvrtdyprime.units=vrt_tendency_units
        self.m_vwndbar_dvrtdyprime=m_vwndbar_dvrtdyprime
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        fileout=fileout.replace(self.filepre,'')
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_vwndbar_dvrtdyprime,fileout)
        if self.archive:
            archive_file(self,fileout)
        #
        ### Calculate m_vwndprime_dvrtdybar
        m_vwndprime_dvrtdybar=-1*vwndprime.data*dvrtdybar.data
        m_vwndprime_dvrtdybar=create_cube(m_vwndprime_dvrtdybar,vrtprime)
        # Attributes
        var_name='m_vwndprime_dvrtdybar'
        long_name=var_name2long_name[var_name]
        m_vwndprime_dvrtdybar.rename(long_name) # not a standard_name
        m_vwndprime_dvrtdybar.var_name=var_name
        m_vwndprime_dvrtdybar.units=vrt_tendency_units
        self.m_vwndprime_dvrtdybar=m_vwndprime_dvrtdybar
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        fileout=fileout.replace(self.filepre,'')
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_vwndprime_dvrtdybar,fileout)
        if self.archive:
            archive_file(self,fileout)
        #
        ### Calculate m_vwndprime_dvrtdyprime
        m_vwndprime_dvrtdyprime=-1*vwndprime*dvrtdyprime
        # Attributes
        var_name='m_vwndprime_dvrtdyprime'
        long_name=var_name2long_name[var_name]
        m_vwndprime_dvrtdyprime.rename(long_name) # not a standard_name
        m_vwndprime_dvrtdyprime.var_name=var_name
        m_vwndprime_dvrtdyprime.units=vrt_tendency_units
        self.m_vwndprime_dvrtdyprime=m_vwndprime_dvrtdyprime
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        fileout=fileout.replace(self.filepre,'')
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_vwndprime_dvrtdyprime,fileout)
        if self.archive:
            archive_file(self,fileout)

    def f_vrt_components(self,level):
        """Calculate and save dvwnddx and m_duwnddy contributions to vorticity.

        dvwndx = vorticity calculated from (0,vwnd)
        m_duwndy = vorticity calculated from (uwnd,0)

        Create attributes:

        dvwnddx
        m_duwnddy
        """
        # Read uwnd and vwnd data for current time block
        self.time1,self.time2=block_times(self,verbose=self.verbose)
        time_constraint=set_time_constraint(self.time1,self.time2,calendar=self.calendar,verbose=self.verbose)
        x1=self.data_in['uwnd_'+str(level)].extract(time_constraint)
        x2=self.data_in['vwnd_'+str(level)].extract(time_constraint)
        self.uwnd=x1.concatenate_cube()
        self.vwnd=x2.concatenate_cube()
        # Create dummy cube of zeroes
        zeros=np.zeros(self.uwnd.data.shape)
        zeros=create_cube(zeros,self.uwnd,new_var_name='dummy')
        # Find value of south2north
        self.south2north=f_south2north(self.uwnd,verbose=self.verbose)
        # Calculate dvwnddx and save
        ww=VectorWind(zeros,self.vwnd)
        self.dvwnddx=ww.vorticity()
        if self.south2north:
            self.dvwnddx=lat_direction(self.dvwnddx,'s2n')
        var_name='dvwnddx'
        long_name=var_name2long_name[var_name]
        self.dvwnddx.rename(long_name) # not a standard_name
        self.dvwnddx.var_name=var_name
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.dvwnddx,fileout)
        if self.archive:
            archive_file(self,fileout)
        # Calculate m_duwnddy and save
        ww=VectorWind(self.uwnd,zeros)
        self.m_duwnddy=ww.vorticity()
        if self.south2north:
            self.m_duwnddy=lat_direction(self.m_duwnddy,'s2n')
        var_name='m_duwnddy'
        long_name=var_name2long_name[var_name]
        self.m_duwnddy.rename(long_name) # not a standard_name
        self.m_duwnddy.var_name=var_name
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.m_duwnddy,fileout)
        if self.archive:
            archive_file(self,fileout)

    def f_div_components(self,level):
        """Calculate and save duwnddx and dvwnddy contributions to divergence.

        duwnddx = divergence calculated from (uwnd,0)
        dvwnddy = divergence calculated from (0,vwnd)

        Create attributes:

        duwnddx
        dvwnddy
        """
        # Read uwnd and vwnd data for current time block
        self.time1,self.time2=block_times(self,verbose=self.verbose)
        time_constraint=set_time_constraint(self.time1,self.time2,calendar=self.calendar,verbose=self.verbose)
        x1=self.data_in['uwnd_'+str(level)].extract(time_constraint)
        x2=self.data_in['vwnd_'+str(level)].extract(time_constraint)
        self.uwnd=x1.concatenate_cube()
        self.vwnd=x2.concatenate_cube()
        # Create dummy cube of zeroes
        zeros=np.zeros(self.uwnd.data.shape)
        zeros=create_cube(zeros,self.uwnd,new_var_name='dummy')
        # Find value of south2north
        self.south2north=f_south2north(self.uwnd,verbose=self.verbose)
        # Calculate duwnddx and save
        ww=VectorWind(self.uwnd,zeros)
        self.duwnddx=ww.divergence()
        if self.south2north:
            self.duwnddx=lat_direction(self.duwnddx,'s2n')
        var_name='duwnddx'
        long_name=var_name2long_name[var_name]
        self.duwnddx.rename(long_name) # not a standard_name
        self.duwnddx.var_name=var_name
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.duwnddx,fileout)
        if self.archive:
            archive_file(self,fileout)
        # Calculate dvwnddy and save
        ww=VectorWind(zeros,self.vwnd)
        self.dvwnddy=ww.divergence()
        if self.south2north:
            self.dvwnddy=lat_direction(self.dvwnddy,'s2n')
        var_name='dvwnddy'
        long_name=var_name2long_name[var_name]
        self.dvwnddy.rename(long_name) # not a standard_name
        self.dvwnddy.var_name=var_name
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.dvwnddy,fileout)
        if self.archive:
            archive_file(self,fileout)

    def f_vrtbudget_combine(self):
        """Combine terms in the vorticity budget.

        Read in individual vorticity budget terms, previously
        calculated using f_vrtbudget().  Combine them together into
        physically useful groups.

        Calculate and create attributes:

        self.vrt_horiz_adv : -u d zeta/dx -v d zeta/dy

        self.vrt_stretch : -zeta D - fD

        self.vrt_tilt : - d omega/dx * dv/dp +d omega/dy * du/dp

        """
        # Read in m_uwnd_dvrtx, m_vwnd_dvrtdy, m_vrt_div, m_ff_div,
        # m_domegadx_dvwnddp, domegady_duwnddp
        # for current time block and assign to attributes
        self.time1,self.time2=block_times(self,verbose=self.verbose)
        time_constraint=set_time_constraint(self.time1,self.time2,calendar=self.calendar,verbose=self.verbose)
        x1=self.data_in['m_uwnd_dvrtdx_'+str(self.level)].extract(time_constraint)
        x2=self.data_in['m_vwnd_dvrtdy_'+str(self.level)].extract(time_constraint)
        x3=self.data_in['m_vrt_div_'+str(self.level)].extract(time_constraint)
        x4=self.data_in['m_ff_div_'+str(self.level)].extract(time_constraint)
        x5=self.data_in['m_domegadx_dvwnddp_'+str(self.level)].extract(time_constraint)
        x6=self.data_in['domegady_duwnddp_'+str(self.level)].extract(time_constraint)
        self.m_uwnd_dvrtdx=x1.concatenate_cube()
        self.m_vwnd_dvrtdy=x2.concatenate_cube()
        self.m_vrt_div=x3.concatenate_cube()
        self.m_ff_div=x4.concatenate_cube()
        self.m_domegadx_dvwnddp=x5.concatenate_cube()
        self.domegady_duwnddp=x6.concatenate_cube()
        # Calculate horizontal advection of vorticity
        self.vrt_horiz_adv=self.m_uwnd_dvrtdx+self.m_vwnd_dvrtdy
        var_name='vrt_horiz_adv'
        long_name=var_name2long_name[var_name]
        self.vrt_horiz_adv.rename(long_name) # not a standard_name
        self.vrt_horiz_adv.var_name=var_name
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.vrt_horiz_adv,fileout)
        if self.archive:
            archive_file(self,fileout)
        # Calculate vortex stretching
        self.vrt_stretch=self.m_vrt_div+self.m_ff_div
        var_name='vrt_stretch'
        long_name=var_name2long_name[var_name]
        self.vrt_stretch.rename(long_name) # not a standard_name
        self.vrt_stretch.var_name=var_name
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.vrt_stretch,fileout)
        if self.archive:
            archive_file(self,fileout)
        # Calculate vortex tilting / twisting
        self.vrt_tilt=self.m_domegadx_dvwnddp+self.domegady_duwnddp
        var_name='vrt_tilt'
        long_name=var_name2long_name[var_name]
        self.vrt_tilt.rename(long_name) # not a standard_name
        self.vrt_tilt.var_name=var_name
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.vrt_tilt,fileout)
        if self.archive:
            archive_file(self,fileout)
        
    def f_overturning_potential(self):
        """Calculate mu at current pressure level.

        Based on Schwendike et al. (2014), equation 1.

        mu is minus the inverse Laplacian of omega.
        """
        # Read data for current time block
        self.time1,self.time2=block_times(self,verbose=self.verbose)
        time_constraint=set_time_constraint(self.time1,self.time2,calendar=self.calendar,verbose=self.verbose)
        x1=self.data_in['omega_'+str(self.level)].extract(time_constraint)
        self.omega=x1.concatenate_cube()
        # Convert omega units to Pa s-1 if not already
        if self.omega.units!='Pa s-1':
            print('Converting omega. Original units {0!s}'.format(self.omega.units))
            self.omega.convert_units('Pa s-1')
            print('Final units {0!s}'.format(self.omega.units))
        # Calculate mu and save
        self.mu=-1*inverse_laplacian(self.omega)
        var_name='mu'
        long_name=var_name2long_name[var_name]
        self.mu.rename(long_name) # not a standard_name
        self.mu.var_name=var_name
        self.mu.units='Pa m2 s-1'
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.mu,fileout)
        if self.archive:
            archive_file(self,fileout)

    def f_overturning_vector_streamfunction(self):
        """Calculate (psioc_lambda,psioc_phi) at current pressure level.

        Based on Schwendike et al. (2014), between equation 1 and 2.

        (psioc_lambda,psioc_phi) is -grad_H(mu).
        """
        # Read data for current time block
        self.time1,self.time2=block_times(self,verbose=self.verbose)
        time_constraint=set_time_constraint(self.time1,self.time2,calendar=self.calendar,verbose=self.verbose)
        x1=self.data_in['mu_'+str(self.level)].extract(time_constraint)
        self.mu=x1.concatenate_cube()
        # Calculate (psioc_lambda,psioc_phi) and save
        x1,x2=gradient(self.mu)
        self.psioc_lambda=-1*x1
        self.psioc_phi=-1*x2
        psioc_units='Pa m s-1' # equivalent to 'kg s-3'
        # psioc_lambda
        var_name='psioc_lambda'
        long_name=var_name2long_name[var_name]
        self.psioc_lambda.rename(long_name) # not a standard_name
        self.psioc_lambda.var_name=var_name
        self.psioc_lambda.units=psioc_units
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.psioc_lambda,fileout)
        if self.archive:
            archive_file(self,fileout)
        # psioc_phi
        var_name='psioc_phi'
        long_name=var_name2long_name[var_name]
        self.psioc_phi.rename(long_name) # not a standard_name
        self.psioc_phi.var_name=var_name
        self.psioc_phi.units=psioc_units
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.psioc_phi,fileout)
        if self.archive:
            archive_file(self,fileout)

    def f_overturning_vector_streamfunction_old(self,method=1):
        """Calculate overturning vector streamfunction at current pressure level.

        Do not use this method. Superseded by f_overturning_vector_streamfunction.

        Based on Schwendike et al. (2014), equation 3.

        psioc_lambda and psioc_phi are components of overturning vector streamfunction.

        psioc_lambda (p) = -\int_{p_s}^p uwndchi(p') dp'
        psioc_phi (p)    = -\int_{p_s}^p vwndchi(p') dp'

        """
        raise ToDoError('Do not use this method. Superseded by f_overturning_vector_streamfunction.')
        # Check levels are correct type and units
        if self.levels.level_type=='plev':
            if self.levels.units=='hPa':
                level_conversion_factor=100 # Convert to Pa (SI units)
            else:
                raise ToDoError('Need to code up for units other than hPa.')
        else:
            raise UserWarning('level_type must be pressure plev: {0!s}'.format(self.levels.level_type))

        # Read data for current time block
        self.time1,self.time2=block_times(self,verbose=self.verbose)
        time_constraint=set_time_constraint(self.time1,self.time2,calendar=self.calendar,verbose=self.verbose)
        x1=self.data_in['uwndchi_'+str(self.level)].extract(time_constraint)
        x2=self.data_in['vwndchi_'+str(self.level)].extract(time_constraint)
        x3=self.data_in['psfc_1'].extract(time_constraint)
        self.uwndchi=x1.concatenate_cube()
        self.vwndchi=x2.concatenate_cube()
        self.psfc=x3.concatenate_cube()
        # Convert psfc units to Pascals if not already
        if self.psfc.units!='Pascals':
            print('Converting psfc. Original units {0!s}'.format(self.psfc.units))
            self.psfc.convert_units('Pascals')
            print('Final units {0!s}'.format(self.psfc.units))
        self.level_width=self.levels.level_widths[self.level]*level_conversion_factor
        self.half_level_below=self.levels.half_levels_below[self.level]*level_conversion_factor
        self.half_level_above=self.levels.half_levels_above[self.level]*level_conversion_factor
        print('level: {0!s}'.format(self.level))
        print('level_width, half_level_below, half_level_above: {0!s}, {1!s}, {2!s}'.format(self.level_width,self.half_level_below,self.half_level_above))
        if method==1:
            # Calculate effective level width dependent on psfc
            # 1) If psfc > half_level_below, use level_width
            # 2) If half_level_below > psfc > half_level_above, 
            #      use psfc-half_level_above
            # 3) If half_level_above > psfc, use 0
            # First, replace missing values of psfc (there are a very small
            #  number in erainterim) with maximum allowed value
            psfc_max_allowed=107000 # Pa
            psfc_capped=np.where(np.greater(self.psfc.data,psfc_max_allowed),psfc_max_allowed,self.psfc.data)
            # Initialise with 2)
            #x1=self.psfc.data-self.half_level_above
            x1=psfc_capped-self.half_level_above
            # Test for 1)
            #x2=np.ma.where(np.ma.greater(self.psfc.data,self.half_level_below),self.level_width,x1)
            x2=np.where(np.greater(psfc_capped,self.half_level_below),self.level_width,x1)
            # Test for 3)
            #x3=np.ma.where(np.ma.greater(self.half_level_above,self.psfc.data,),0,x2)
            x3=np.where(np.greater(self.half_level_above,psfc_capped,),0,x2)
            #print('level width min, max: {0!s}, {1!s}'.format(x3.data.min(),x3.data.max()))
            print('level width min, max: {0!s}, {1!s}'.format(x3.min(),x3.max()))
            self.effective_level_width=x3
        elif method==2:
            # Do not use actual psfc
            # Instead, integrate from a constant 1013 hPa surface pressure.
            self.effective_level_width=self.level_width
        else:
            raise UserWarning('Invalid value of method: {0!s}'.format(method))
        # Calculate increments to psioc_lambda and psioc_phi at this level
        # Note, although there is a minus sign in the equation for the
        #   line below (integrating equation 3 of Schwendike et al., 2014), 
        #   we are integrating from the surface upwards, i.e., in the 
        #   direction of negative p, hence we need another minus sign. Result 
        #   is there is no minus sign in line of code below.
        self.psioc_lambda_inc=self.uwndchi*self.effective_level_width
        self.psioc_phi_inc=self.vwndchi*self.effective_level_width
        psioc_units='Pa m s-1' # equivalent to 'kg s-3'
        self.psioc_lambda_inc.units=psioc_units
        self.psioc_phi_inc.units=psioc_units
        # Add increments to psioc_lambda and psioc_phi from previous level
        if self.level==self.levels.levels[0]:
            print('First level. psioc = increment.')
            self.psioc_lambda=self.psioc_lambda_inc
            self.psioc_phi=self.psioc_phi_inc
        else:
            print('Not first level. Increment added to psioc.')
            self.psioc_lambda=self.psioc_lambda+self.psioc_lambda_inc
            self.psioc_phi=self.psioc_phi+self.psioc_phi_inc
        # Save psioc_lambda
        self.file_data_out=os.path.join(self.basedir,self.source,'std','VAR_NAME_'+str(self.level)+self.filepre+'_'+self.wildcard+'.nc')
        var_name='psioc_lambda'
        long_name=var_name2long_name[var_name]
        self.psioc_lambda.rename(long_name) # not a standard_name
        self.psioc_lambda.var_name=var_name
        self.psioc_lambda.units=psioc_units
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.psioc_lambda,fileout)
        if self.archive:
            archive_file(self,fileout)
        # Save psioc_phi
        var_name='psioc_phi'
        long_name=var_name2long_name[var_name]
        self.psioc_phi.rename(long_name) # not a standard_name
        self.psioc_phi.var_name=var_name
        self.psioc_phi.units=psioc_units
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.psioc_phi,fileout)
        if self.archive:
            archive_file(self,fileout)

    def f_omega_decomposition(self):
        """Calculate omega_lambda and omega_phi at current pressure level.

        Based on Schwendike et al. (2014), equations 4 and 5.

        psioc_lambda and psioc_phi are components of overturning vector streamfunction.

        omega_lambda and omega_phi is decomposition of omega into 
        eastward and northward overturning circulations

        omega_lambda = (1 / a cos phi) d psi_lambda /d lambda
                     = "x" component of grad.(psi_lambda,0)

        omega_phi = (1 / a cos phi) d / d phi (cos phi psi_phi)
                     = "y" component of grad.(0,psi_phi)
        """
        # Read data for current time block
        self.time1,self.time2=block_times(self,verbose=self.verbose)
        time_constraint=set_time_constraint(self.time1,self.time2,calendar=self.calendar,verbose=self.verbose)
        x1=self.data_in['psioc_lambda_'+str(self.level)].extract(time_constraint)
        x2=self.data_in['psioc_phi_'+str(self.level)].extract(time_constraint)
        self.psioc_lambda=x1.concatenate_cube()
        self.psioc_phi=x2.concatenate_cube()
        # Create dummy cube of zeroes
        zeros=np.zeros(self.psioc_lambda.data.shape)
        zeros=create_cube(zeros,self.psioc_lambda,new_var_name='dummy')
        # Calculate omega_lambda and save
        ww=VectorWind(self.psioc_lambda,zeros)
        self.omega_lambda=ww.divergence()
        var_name='omega_lambda'
        long_name=var_name2long_name[var_name]
        self.omega_lambda.rename(long_name) # not a standard_name
        self.omega_lambda.var_name=var_name
        omega_units='Pa s-1'
        self.omega_lambda.units=omega_units
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.omega_lambda,fileout)
        if self.archive:
            archive_file(self,fileout)
        # Calculate omega_phi and save
        ww=VectorWind(zeros,self.psioc_phi)
        self.omega_phi=ww.divergence()
        var_name='omega_phi'
        long_name=var_name2long_name[var_name]
        self.omega_phi.rename(long_name) # not a standard_name
        self.omega_phi.var_name=var_name
        self.omega_phi.units=omega_units
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.omega_phi,fileout)
        if self.archive:
            archive_file(self,fileout)

    def f_mass_flux_decomposition(self):
        """Calculate mf_lambda and mf_phi at current pressure level.

        Based on Schwendike et al. (2014), equation 6.

        omega_lambda and omega_phi is decomposition of omega into 
        eastward and northward overturning circulations

        mf_lambda and mf_phi is decomposition of upward mass flux into 
        eastward and northward overturning circulations

        mg_lambda = -omega_lambda cos phi / g
        mg_phi    = -omega_phi    cos phi / g

        """
        # Read data for current time block
        self.time1,self.time2=block_times(self,verbose=self.verbose)
        time_constraint=set_time_constraint(self.time1,self.time2,calendar=self.calendar,verbose=self.verbose)
        x1=self.data_in['omega_lambda_'+str(self.level)].extract(time_constraint)
        x2=self.data_in['omega_phi_'+str(self.level)].extract(time_constraint)
        self.omega_lambda=x1.concatenate_cube()
        self.omega_phi=x2.concatenate_cube()
        # Calculate cos phi array
        lats=self.omega_lambda.coord('latitude').points
        nlats=len(lats)
        cosphi=np.cos(lats*np.pi/180)
        shape=self.omega_lambda.shape
        lat_index=shape.index(nlats)
        shape=list(shape)
        for index in range(len(shape)):
            if index!=lat_index:
                shape[index]=1
        shape=tuple(shape)
        cosphi=cosphi.reshape(shape)
        # Calculate mf_lambda and save
        self.mf_lambda=-1*self.omega_lambda*cosphi/info.gg
        var_name='mf_lambda'
        long_name=var_name2long_name[var_name]
        self.mf_lambda.rename(long_name) # not a standard_name
        self.mf_lambda.var_name=var_name
        mf_units='kg m-2 s-1'
        self.mf_lambda.units=mf_units
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.mf_lambda,fileout)
        if self.archive:
            archive_file(self,fileout)
        # Calculate mf_phi and save
        self.mf_phi=-1*self.omega_phi*cosphi/info.gg
        var_name='mf_phi'
        long_name=var_name2long_name[var_name]
        self.mf_phi.rename(long_name) # not a standard_name
        self.mf_phi.var_name=var_name
        self.mf_phi.units=mf_units
        fileout=self.file_data_out.replace('VAR_NAME',var_name)
        fileout=replace_wildcard_with_time(self,fileout)
        print('fileout: {0!s}'.format(fileout))
        iris.save(self.mf_phi,fileout)
        if self.archive:
            archive_file(self,fileout)

    def f_erainterim_calc_levels(self):
        """Calculate ERA-Interim pressure half and full levels from surface pressure.

        Based on Berrisford et al. (2011)

        """
        # Read psfc for current time block
        self.time1,self.time2=block_times(self,verbose=self.verbose)
        time_constraint=set_time_constraint(self.time1,self.time2,calendar=self.calendar,verbose=self.verbose)
        x1=self.data_in['psfc_'+str(self.level)].extract(time_constraint)
        self.psfc=x1.concatenate_cube()
        # Calculate pressure of half and full levels for model levels
        nlevel=info.nlevel_erainterim_model
        print('nlevel: {0!s}'.format(nlevel))
        # Level 1
        ii=1
        akmhalf=info.erainterim_abkmalf[ii][0]
        bkmhalf=info.erainterim_abkmalf[ii][1]
        print('ii,akmhalf,bkmhalf: {0!s}, {1!s}, {2!s}'.format(ii,akmhalf,bkmhalf))
        xx=akmhalf+bkmhalf*self.psfc
        xxmin=xx.data.min()
        xxmean=xx.data.mean()
        xxmax=xx.data.max()
        print('xxmin,xxmean,xxmax: {0!s}, {1!s}, {2!s}'.format(xxmin,xxmean,xxmax))
        # Level 2
        xxm1=xx
        ii=2
        akmhalf=info.erainterim_abkmalf[ii][0]
        bkmhalf=info.erainterim_abkmalf[ii][1]
        print('ii,akmhalf,bkmhalf: {0!s}, {1!s}, {2!s}'.format(ii,akmhalf,bkmhalf))
        xx=akmhalf+bkmhalf*self.psfc
        xxmin=xx.data.min()
        xxmean=xx.data.mean()
        xxmax=xx.data.max()
        print('xxmin,xxmean,xxmax: {0!s}, {1!s}, {2!s}'.format(xxmin,xxmean,xxmax))
        pdb.set_trace()
        #for ii in range(2,nlevel+1):
        for ii in range(2,3+1):
        #for ii in range(51,52+1):
            akmhalf=info.erainterim_abkmalf[ii][0]
            bkmhalf=info.erainterim_abkmalf[ii][1]
            print('ii,akmhalf,bkmhalf: {0!s}, {1!s}, {2!s}'.format(ii,akmhalf,bkmhalf))
            xx=akmhalf+bkmhalf*self.psfc
            xxmin=xx.data.min()
            xxmean=xx.data.mean()
            xxmax=xx.data.max()
            print('xxmin,xxmean,xxmax: {0!s}, {1!s}, {2!s}'.format(xxmin,xxmean,xxmax))


#==========================================================================

class Planet(object):

    """Planetary parameters.

    Default is Planet Earth.

    Attributes:

    Specified base attributes

    self.arad : Planet radius (m)
    self.gg : Gravitational acceleration at surface (m s-2)
    self.omega :  Sidereal rotation rate (s-1)
    self.gascon : Specific gas constant of atmosphere (J kg-1 K-1)
    self.cp : Specific heat capacity of atmosphere at constant pressure (J kg-1 K-1)

    Derived attributes

    self.f0 : Coriolis parameter at 45 degrees lat (s-1)
    self.beta0 : beta at equator (m-1 s-1)
    self.circum : planet circumference (m)
    """

    def __init__(self,title='Earth',verbose=True):
        """Initialise."""
        self.title=title
        self.verbose=verbose
        if self.title=='Earth':
            self.arad=6370000 # m
            self.gg=9.81 # m s-2
            self.omega=7.292e-5 # s-1
            self.gascon=287 # J kg-1 K-1
            self.cp=1004 # J kg-1 K-1
        else:
            raise ToDoError('Code up for other planets.')
        self.f0=np.sqrt(2)*self.omega
        self.beta0=2*self.omega/self.arad
        self.circum=2*np.pi*self.arad
        if self.verbose:
            print(self)
            
    def __repr__(self):
        return 'Planet(title={0.title!r}, verbose={0.verbose!r})'.format(self)
    
    def __str__(self):
        if self.verbose:
            ss=h1a+'Planet instance \n'+\
                'title: {0.title!s} \n'+\
                'arad: {0.arad!s} m \n'+\
                'gg: {0.gg!s} m s-2 \n'+\
                'omega: {0.omega!s} s-1 \n'+\
                'gascon: {0.gascon!s} J kg-1 K-1 \n'+\
                'cp: {0.cp!s} J kg-1 K-1 \n'+\
                'f0: {0.f0!s} s-1 \n'+\
                'beta0: {0.beta0!s} m-1 s-1 \n'+\
                'circum: {0.circum!s} m \n'+h1b
            return ss.format(self)
        else:
            return self.__repr__()    

#==========================================================================

class XYGrid(object):

    """Regularly spaced 2-d rectangular grid.
    
    Basic grid is initially an empty (zeroes) 2-d numpy array.  Populate
    the grid with numbers externally.

    Two alternative sets of axes are used: (x,y) and (lon,lat).  All
    calculations on the grid are done using the (x,y) axes.  The
    (lon,lat) axes are only for presentational purposes.

    Attributes:
    
    self.nx=self.nlon, self.ny=self.nlat : integer number of grid
    points in x and y directions.
    
    self.deltax, self.deltay : float spacing for x and y grid points
    
    self.deltalon, self.deltalat : equivalent float spacing for lon
    and lat grid points
    
    self.x1, self.y1 : float values of first points on x and y axes
    
    self.lon1, self.lat1 : equivalent float values of first points on
    lon and lat axes

    self.xaxis, self.yaxis : numpy arrays of x and y values
    
    self.lonaxis, self.lataxis : equivalent numpy arrays of lon and
    lat values
    
    self.xperiodic, self.yperiodic : Boolean flags for whether x and y
    axes are periodic (wrapping around).
    
    """
    
    def __init__(self,**descriptor):
        """Initialise from descriptor dictionary.

        Compulsory keywords: 'verbose','var_name','nx','ny','type_init',
        'units','xperiodic','yperiodic'.

        Optional keywords: ['deltalon','deltalat','lon1','lat1'] or
        ['deltax','deltay','x1','y1'].
        """
        self.__dict__.update(descriptor)
        self.nlon=self.nx
        self.nlat=self.ny
        if self.type_init not in ['lonlat','xy']:
            raise UserWarning('Invalid initial type.')
        # Create initial 2-d grid of zeros in numpy array
        self.grid=np.zeros((self.ny,self.nx))
        self.planet=Planet()
        if self.type_init=='lonlat':
            # Create x-y equivalents
            self.deltax=(self.deltalon/360)*self.planet.circum
            self.deltay=(self.deltalat/360)*self.planet.circum
            self.x1=(self.lon1/360)*self.planet.circum
            self.y1=(self.lat1/360)*self.planet.circum
        elif self.type_init=='xy':
            # Create lon-lat equivalents
            self.deltalon=(self.deltax/self.planet.circum)*360
            self.deltalat=(self.deltay/self.planet.circum)*360
            self.lon1=(self.x1/self.planet.circum)*360
            self.lat1=(self.y1/self.planet.circum)*360
        # Create lon-axis
        self.lon2=self.lon1+(self.nlon-1)*self.deltalon
        epsilon=self.deltalon/100
        self.lonaxis=np.arange(self.lon1,self.lon2+epsilon,self.deltalon)
        # Create lat-axis
        self.lat2=self.lat1+(self.nlat-1)*self.deltalat
        self.lataxis=np.arange(self.lat1,self.lat2+epsilon,self.deltalat)
        # Create x-axis
        self.x2=self.x1+(self.nx-1)*self.deltax
        epsilon=self.deltax/100
        self.xaxis=np.arange(self.x1,self.x2+epsilon,self.deltax)
        # Create y-axis
        self.y2=self.y1+(self.ny-1)*self.deltay
        self.yaxis=np.arange(self.y1,self.y2+epsilon,self.deltay)
        if self.verbose:
            print(self)
            
    def __repr__(self):
        return 'XYGrid({0.type_init!r},verbose={0.verbose!r})'.format(self)
    
    def __str__(self):
        if self.verbose:
            ss=h1a+'XYGrid instance \n'+\
                'type_init: {0.type_init!s} \n'+\
                'nlon,nlat: {0.nlon!s},{0.nlat!s} \n'+\
                'deltalon,deltalat: {0.deltalon!s},{0.deltalat!s} \n'+\
                'lon1,lon2: {0.lon1!s},{0.lon2!s} \n'+\
                'lat1,lat2: {0.lat1!s},{0.lat2!s} \n'+\
                'nx,ny: {0.nx!s},{0.ny!s} \n'+\
                'deltax,deltay: {0.deltax!s},{0.deltay!s} \n'+\
                'x1,x2: {0.x1!s},{0.x2!s} \n'+\
                'y1,y2: {0.y1!s},{0.y2!s} \n'+\
                'xperiodic,yperiodic: {0.xperiodic!s},{0.yperiodic!s} \n'+h1b
            return ss.format(self)
        else:
            return self.__repr__()

    def x_gradient(self):
        """Calculate x-gradient of self.grid.

        Creates attribute self.x_gradient."""

    def y_gradient(self):
        """Calculate y-gradient of self.grid.

        Creates attribute self.y_gradient."""

    def create_xy_cube(self):
        """Create iris cube from x-y grid.

        Creates attribute self.xy_cube."""

    def create_lonlat_cube(self):
        """Create iris cube from lon-lat grid.

        Creates attribute self.lonlat_cube."""
        loncoord=iris.coords.DimCoord(self.lonaxis,standard_name='longitude',units='degree_east',circular=self.xperiodic)
        latcoord=iris.coords.DimCoord(self.lataxis,standard_name='latitude',units='degree_north')
        dim_coords=[[latcoord,0],[loncoord,1]]
        newcube=iris.cube.Cube(self.grid,units=self.units,dim_coords_and_dims=dim_coords)
        long_name=var_name2long_name[self.var_name]
        newcube.rename(long_name)
        newcube.var_name=self.var_name
        self.lonlatcube=newcube
        
#==========================================================================

class WheelerKiladis(object):

    """Wheeler-Kiladis diagnostics of 2-D longitude-time data..

    Called from wheelerkiladis.py.

    Selected attributes:

    self.data_hov : iris cube of input 2-D longitude-time data

    self.data_hovfftWK : iris cube of 2-D FFT of self.hov (complex
    data, so cannot be saved as netcdf file)

    self.data_hovfftWKabs : iris cube of magnitude of 2-D FFT of
    self.hov

    self.data_hovWKfilt : inverse 2-D FFT of self.hovfft, i.e., 2-D
    longitude-time Hovmoller of wave-filtered data

    """

    def __init__(self,**descriptor):
        """Initialise from descriptor dictionary.

        Compulsory keywords: 'verbose','source','var_name','level',
        'basedir','filepre', 'time1', 'time2', 'lat1', 'lat2'.
        """
        self.__dict__.update(descriptor)
        self.descriptor=descriptor
        self.name=var_name2long_name[self.var_name]
        source_info(self)
        # Input data
        self.filein1=os.path.join(self.basedir,self.source,'processed',self.var_name+'_'+str(self.level)+self.filepre+'_ss_lat_'+str(self.lat1)+'_'+str(self.lat2)+'_'+self.wildcard+'.nc')
        self.data_in=iris.load(self.filein1,self.name)
        # Output files
        ss='_lat_'+str(self.lat1)+'_'+str(self.lat2)+'_'+str(self.time1)[:10]+'_'+str(self.time2)[:10]
        self.file_hov=os.path.join(self.basedir,self.source,'processed',self.var_name+'_'+str(self.level)+self.filepre+'_hov'+ss+'.nc')
        self.file_hovfftWKabs=os.path.join(self.basedir,self.source,'processed',self.var_name+'_'+str(self.level)+self.filepre+'_hovfftWKabs'+self.wave_type+ss+'.nc')
        self.file_hovWKfilt=os.path.join(self.basedir,self.source,'processed',self.var_name+'_'+str(self.level)+self.filepre+'_hovWKfilt'+self.wave_type+ss+'.nc')
        if self.verbose:
            print(self)
        
    def __repr__(self):
        return 'WheelerKiladis({0.descriptor!r},verbose={0.verbose!r})'.format(self)

    def __str__(self):
        if self.verbose==2:
            ss=h1a+'WheelerKiladis instance \n'+\
                'var_name: {0.var_name!s} \n'+\
                'level: {0.level!s} \n'+\
                'source: {0.source!s} \n'+\
                'data_in: {0.data_in!s} \n'+\
                'filein1: {0.filein1!s} \n'+\
                'file_hov: {0.file_hov!s} \n'\
                'file_hovWKfilt: {0.file_hovWKfilt!s} \n'+h1b
            return ss.format(self)
        else:
            return 'WheelerKiladis of '+self.source+' '+self.var_name+str(self.level)

    def f_hovmoller(self,mode='calculate'):
        """Create and write single cube 2-D Hovmoller, or read pre-existing.

        Inputs:

        If <mode> is 'calculate', create a single cube 2-D Hovmoller
        self.data_hov from the input data between self.time1 and
        self.time2, and save to file.

        If <mode> is 'read', read in pre-existing single cube 2-D
        Hovmoller self.data_hov from file.

        If <mode> is 'analytical', create an analytical test Hovmoller
        with single propagating wave.

        self.<hovmoller_params> : dictionary of parameters used to
        create the Hovmoller, related to preparing it for future 2-D
        FFT.

        Output:

        Create attribute self.<data_hov> : iris cube of Hovmoller data

        """
        if mode=='calculate':
            print('Creating single cube 2-D Hovmoller.')
            if self.calendar=='gregorian':
                pass
            elif self.calendar=='360_day':
                raise ToDoError('Need to think how this will work with 360_day calendar.')
            else:
                raise UserWarning('Invalid calendar.')
            time_constraint=set_time_constraint(self.time1,self.time2,calendar=self.calendar,verbose=self.verbose)
            xx1=self.data_in.extract(time_constraint)
            self.data_hov=xx1.concatenate_cube()
            time_coord=self.data_hov.coord('time')
            print(time_coord[0])
            print(time_coord[-1])
            print('Saving single cube 2-D Hovmoller')
            iris.save(self.data_hov,self.file_hov)
            if self.archive:
                archive_file(self,self.file_hov)
        elif mode=='read':
            print('Reading single cube 2-D Hovmoller')
            self.data_hov=iris.load_cube(self.file_hov)
        elif mode=='analytical':
            print('Creating analytical test Hovmoller')
            # This requires observed Hovmoller to have been previously
            # calculated and saved. That is read in here, and overwritten
            # with a single propagating wave to test FFT and rearrangement
            # code later
            xx1=iris.load_cube(self.file_hov)
            self.nt=xx1.coord('time').points.shape[0]
            self.nx=xx1.coord('longitude').points.shape[0]
            self.ntd2=int(self.nt/2)
            self.nxd2=int(self.nx/2)
            deltat=1
            deltax=1
            amp=1
            iharm=36 # frequency (harmonic) of input wave(>=0)
            ss=3 # zonal wavenumber of input wave (pos or neg)
            if not(0<=iharm<=self.ntd2):
                raise UserWarning('iharm out of range.')
            if not(-self.nxd2<=ss<=self.nxd2):
                raise UserWarning('ss out of range.')
            if iharm==0:
                omega=0
            else:
                omega=2*np.pi/(self.nt*deltat/iharm)
            if ss==0:
                kk=0
            else:
                kk=2*np.pi/(self.nx*deltax/ss)
            phi=0.11*2*np.pi
            tt=np.arange(self.nt).reshape((self.nt,1))
            ones=np.ones((1,self.nx))
            tt=np.dot(tt,ones)
            ones=np.ones((self.nt,1))
            xx=np.arange(self.nx).reshape((1,self.nx))
            xx=np.dot(ones,xx)
            hov=amp*np.cos(kk*xx-omega*tt-phi)
            self.data_hov=create_cube(hov,xx1)
        else:
            raise UserWarning('Invalid mode.')
        self.nt=self.data_hov.coord('time').points.shape[0]
        self.nx=self.data_hov.coord('longitude').points.shape[0]
        self.ntd2=int(self.nt/2)
        self.nxd2=int(self.nx/2)
        print('data_hov: {0.data_hov!s}'.format(self))
        if self.hovmoller_params['cosine_tapering']:
            raise ToDoError('The cosine tapering code needs careful checking, and is probably not necessary anyway.')
            frac=self.hovmoller_params['cosine_tapering_fraction']
            print('Applying cosine tapering in time to first and last {0!s} fraction of Hovmoller.'.format(frac))
            taper=np.ones((self.nt))
            nn=len(taper)
            ntaper=int(frac*nn)
            nper=2*ntaper
            for ii in range(ntaper):
                taper[ii]=1-(0.5*(np.cos(2*np.pi*ii/nper)+1))
                taper[nn-1-ii]=1-(0.5*(np.cos(2*np.pi*ii/nper)+1))
            ones=np.ones((1,self.nx))
            taper=taper.reshape((nn,1))
            taper=np.dot(taper,ones)
            print('taper: {0!s}'.format(taper))
            self.taper=taper
            xx=self.data_hov.data.mean()
            print('Mean of Hovmoller before tapering: {0!s}'.format(xx))
            self.data_hov=create_cube(taper*self.data_hov.data,self.data_hov)
            cm=iris.coords.CellMethod('point',coords=['time'],comments='cosine tapered in time over first and last '+str(frac)+' fraction of data')
            self.data_hov.add_cell_method(cm)
            xx=self.data_hov.data.mean()
            print('Mean of Hovmoller after tapering: {0!s}'.format(xx))
            print('data_hov: {0.data_hov!s}'.format(self))

    def f_fft(self):
        """Calculate 2-D FFT of 2-D Hovmoller.

        Create a 2-D FFT iris cube self.data_hovfftWK from 2-D
        Hovmoller self.data_hov. Note this is an array of complex
        numbers which cannot be saved to a netcdf file (but could save
        real and imaginary parts separately if needed). So just save
        the absolute values (magnitudes of the complex numbers) in
        self.data_hovfftWKabs to a netcdf file for reference, but
        retain the full complex number array for future work.

        Notes on ouput from np.fft.fft2 and rearrangement of this
        output here to be consistent with Wheeler-Kiladis diagrams.

        Input to fft2 is (nt,nx) Hovmoller. We need both nt and nx to
        be even for subsequent rearrangement.

        Output from fft2 is a 2-D (nt,nx) array.  Arrangment of
        (complex) coefficients is shown by the following. Because the
        input to fft2 is real, there is a symmetry in some of the
        outputs. We could use rfft2 which only outputs approximately
        half the number of coefficients to avoid duplication. However,
        the output of rfft2 is complex, so the filtered version of
        this complex FFT cannot be fed back into irfft2 to compute the
        inverse FFT as irfft2 requires real input. So have to use the
        full fft2.

        For example with nt=6 and nx=8 (and a Nyquist frequency of
        harmonic 3, and Nyquist wavenumber of 4), the output of fft2
        is:

        In the following table, each entry is (iharm,ss) where iharm
        is iharm'th time harmonic and ss is the ss'th zonal
        wavenumber. In this representation, there are only positive
        time harmonics, but ss can be positive or negative
        (corresponding to eastward and westward propagating waves).

        NB The values in the upper left quadrant are the complex
        conjugates of the corresponding coefficients in the lower
        right quadrant (i.e., not actually identical). Note that the
        imaginary parts of the coefficients correspond to a sine wave
        test case sin(kx-omega t), and the real parts correspond to a
        cosine wave test case cos(kx -omega t), so make sure to test
        with both!

        Similarly the values in the upper right quadrant are the
        complex conjugates of the corresponding coefficients in the
        lower left quadrant.

        index
        5,nt-1   | 1,0  1,1  1,2  1,3  1,4 1,-3 1,-2, 1,-1
          4      | 2,0  2,1  2,2  2,3  2,4 2,-3 2,-2  2,-1
        3,nt/2   | 3,0  3,1  3,2  3,3  3,4 3,3  3,2   3,1
          2      | 2,0  2,-1 2,-2 2,-3 2,4 2,3  2,2   2,1
          1      | 1,0  1,-1 1,-2 1,-3 1,4 1,3  1,2   1,1
          0      | 0,0  0,1  0,2  0,3  0,4 0,3  0,2   0,1
                 |________________________________________
                    0    1    2    3    4   5    6     7  x index
                    0                  nx/2           nx-1

        Note that e.g., 1,4 and 1,-4 are identical, as ss=4 is Nyquist
        wavenumber. Also, e.g., 3,2 and 3,-2 are identical as iharm=3
        is Nyquist harmonic.
              
        We only need the left-hand half of this array, i.e., all the t
        indices, and the x indices up to the Nyquist wavenumber (index
        4) in this case. The right hand half is a duplication.

        This needs to be rearranged into a form consistent with the
        Wheeler-Kiladis diagram, a (nt/2+1,nx+1) array:

        iharm  index  | 
        3,nt/2 3,nt/2 | 3,4 3,3  3,2  3,1  3,0 3,1 3,2 3,3 3,4
          2     2     | 2,4 2,-3 2,-2 2,-1 2,0 2,1 2,2 2,3 2,4
          1     1     | 1,4 1,-3 1,-2 1,-1 1,0 1,1 1,2 1,3 1,4
          0     0     | 0,4 0,3  0,2  0,1  0,0 0,1 0,2 0,3 0,4
        --------------|________________________________________
             ss          -4 -3   -2   -1    0   1   2   3   4 
                       -nx/2     -2   -1    0   1   2      nx/2
           index         0   1    2    3    4   5   6   7   8
                         0                 nx/2            nx

        Output. Creates attributes:

        self.data_hovfft : wavenumber-frequency 2-D complex numpy
        array output from np.fft.fft2. Shape (nt,nx)

        self.data_hovfftWK : wavenumber-frequency 2-D complex numpy
        array, in Wheeler-Kiladis format. Shape (nt/2 +1,nx+1)

        self.data_hovfftWKabs : magnitude of data_hovfftWK. Hence
        real. Stored as iris cube.

        """
        print('Creating 2-D FFT of 2-D Hovmoller.')
        # Check that self.data_hov is a (t,x) cube
        if self.data_hov.coords()[0].name()!='time':
            raise UserWarning('Dimension 0 must be time.')
        if self.data_hov.coords()[1].name()!='longitude':
            raise UserWarning('Dimension 0 must be longitude.')
        # Check there are an even number of times and longitudes
        tcoord=self.data_hov.coord('time')
        if divmod(self.nt,2)[1]==0:
            print('Ok, nt is even: {0.nt!s}'.format(self))
        else:
            raise UserWarning('Must be an even number of times: {0.nt!s}'.format(self))
        xcoord=self.data_hov.coord('longitude')
        if divmod(self.nx,2)[1]==0:
            print('Ok, nx is even: {0.nx!s}'.format(self))
        else:
            raise UserWarning('Must be an even number of longitudes: {0.nx!s}'.format(self))
        # Calculate 2D FFT
        # NB At first inspection, would make sense to use the numpy.fft.rfft2
        # method as the input data (Hovmoller) is real.
        # However, in practice we cannot do this.
        # This is because the output FFT is complex. We
        # need to filter this then transform it back using an inverse FFT.
        # The complementary irfft2 method also requires real data so cannot
        # be used.
        # Hence we need to use the general fft2 and ifft2 methods.
        #self.data_hovfft=np.fft.fft2(self.data_hov.data,norm='ortho')
        self.data_hovfft=np.fft.fft2(self.data_hov.data)
        print('Output from ff2 data_hovfft.shape: {0!s}'.format(self.data_hovfft.shape))
        # Rearrange 2D FFT into form consistent with Wheeler-Kiladis diagrams
        x2=self.fft2_numpy_to_WK(self.data_hovfft)
        # Create frequency (and alternative harmonic and angular frequency) 
        # coordinate
        freq_val=np.fft.fftfreq(self.nt)
        freq_val=freq_val[:self.ntd2+1]
        freq_val[-1]=0.5 # Switch from -0.5 to 0.5
        harm_val=np.arange(x2.shape[0])
        harm_units='1'
        freq_units='cycles per '+self.frequency
        if self.frequency=='d':
            factor=1/86400
        elif self.frequency=='6h':
            factor=(24/6)/86400
        elif self.frequency=='3h':
            factor=(24/3)/86400
        else:
            raise ToDoError('Generalise code for other source frequency.')
        omegaf_val=2*np.pi*freq_val*factor
        omegaf_units='rad s-1'
        freq_coord=iris.coords.DimCoord(freq_val,var_name='freq',long_name=var_name2long_name['freq'],units=freq_units)
        harm_coord=iris.coords.DimCoord(harm_val,var_name='harm',long_name=var_name2long_name['harm'],units=harm_units)
        omegaf_coord=iris.coords.DimCoord(omegaf_val,var_name='omegaf',long_name=var_name2long_name['omegaf'],units=omegaf_units)
        # Create zonal wavenumber (and alternative kk wavenumber m-1) coordinate
        ss_val=np.arange(-self.nxd2,self.nxd2+1)
        ss_units='1'
        ss_coord=iris.coords.DimCoord(ss_val,var_name='ss',long_name=var_name2long_name['ss'],units=ss_units)
        planet=Planet()
        kk_val=ss_val/planet.arad
        kk_units='m-1'
        kk_coord=iris.coords.DimCoord(kk_val,var_name='kk',long_name=var_name2long_name['kk'],units=kk_units)
        # Decide on dim and aux coords
        dim_coords_and_dims=[[harm_coord,0],[ss_coord,1]]
        aux_coords_and_dims=[[freq_coord,0],[omegaf_coord,0],[kk_coord,1]]
        # Create iris cube of rearranged 2-D FFT 
        # (array of complex numbers)
        x3=iris.cube.Cube(x2,units='',attributes=self.data_hov.attributes,cell_methods=self.data_hov.cell_methods,dim_coords_and_dims=dim_coords_and_dims,aux_coords_and_dims=aux_coords_and_dims)
        x3.var_name=self.data_hov.var_name+'_fft'
        self.data_hovfftWK=x3
        print('data_hovfftWK: {0.data_hovfftWK!s}'.format(self))
        # Create iris cube of magnitude rearranged 2-D FFT 
        # (array of real numbers which can be saved as netcdf file)
        x3=iris.cube.Cube(np.abs(x2),units='',attributes=self.data_hov.attributes,cell_methods=self.data_hov.cell_methods,dim_coords_and_dims=dim_coords_and_dims,aux_coords_and_dims=aux_coords_and_dims)
        x3.var_name=self.data_hov.var_name+'_fftabs'
        self.data_hovfftWKabs=x3
        iris.save(self.data_hovfftWKabs,self.file_hovfftWKabs)
        if self.archive:
            archive_file(self,self.file_hovfftWKabs)
        if self.verbose==2:
            print('write_cube: {0.file_hovfftWKabs!s}'.format(self))

    def fft2_numpy_to_WK(self,x1):
        """Rearrange output of fft2 into Wheeler-Kiladis format.

        Inputs:

        <x1> : numpy array, output from numpy.fft.fft2. Shape(nt,nx)

        Output:

        <x2> : numpy array, 'rearranged' version of <x1>, not same
        size. Shape (nt/2 +1,nx+1).

        """
        # Create empty array x2 for Wheeler-Kiladis arrangment, then overwrite
        x2=1e20*np.ones((self.ntd2+1,self.nx+1),dtype=np.cdouble)
        # Convention: it,ix in following code refers to index in *source* array
        # Copy to ss<= 0 coefficients in x2
        """
        index
        5,nt-1   |                                          
          4      |                                        
        3,nt/2   | 3,0  3,1  3,2  3,3  3,4               
          2      | 2,0  2,-1 2,-2 2,-3 2,4                 
          1      | 1,0  1,-1 1,-2 1,-3 1,4                 
          0      | 0,0  0,1  0,2  0,3  0,4               
                 |________________________________________
                    0    1    2    3    4   5    6     7  x index
                    0                  nx/2           nx-1

                             copy to

        iharm  index  | 
        3,nt/2 3,nt/2 | 3,4 3,3  3,2  3,1  3,0                
          2     2     | 2,4 2,-3 2,-2 2,-1 2,0                
          1     1     | 1,4 1,-3 1,-2 1,-1 1,0                  
          0     0     | 0,4 0,3  0,2  0,1  0,0                 
        --------------|________________________________________
             ss          -4 -3   -2   -1    0   1   2   3   4 
                       -nx/2     -2   -1    0   1   2      nx/2
           index         0   1    2    3    4   5   6   7   8
                         0                 nx/2            nx"""
        for it in range(self.ntd2+1):
            for ix in range(self.nxd2+1):
                x2[it,self.nxd2-ix]=x1[it,ix]
        # Copy to iharm=0, ss>0 coefficients in x2
        """
        index
        5,nt-1   |                                        
          4      |                                        
        3,nt/2   |                                        
          2      |                                         
          1      |                                       
          0      |      0,1  0,2  0,3  0,4               
                 |________________________________________
                    0    1    2    3    4   5    6     7  x index
                    0                  nx/2           nx-1

                             copy to

        iharm  index  | 
        3,nt/2 3,nt/2 |                                       
          2     2     |                                          
          1     1     |                                       
          0     0     |                        0,1 0,2 0,3 0,4
        --------------|________________________________________
             ss          -4 -3   -2   -1    0   1   2   3   4 
                       -nx/2     -2   -1    0   1   2      nx/2
           index         0   1    2    3    4   5   6   7   8
                         0                 nx/2            nx"""
        it=0
        for ix in range(1,self.nxd2+1):
            x2[it,self.nxd2+ix]=x1[it,ix]
        # Copy to iharm>0, ss>0 coefficients in x2
        """
        index
        5,nt-1   |      1,1  1,2  1,3  1,4                
          4      |      2,1  2,2  2,3  2,4                
        3,nt/2   |      3,1  3,2  3,3  3,4               
          2      |                                       
          1      |                                        
          0      |                                        
                 |________________________________________
                    0    1    2    3    4   5    6     7  x index
                    0                  nx/2           nx-1

                             copy to

        iharm  index  | 
        3,nt/2 3,nt/2 |                        3,1 3,2 3,3 3,4
          2     2     |                        2,1 2,2 2,3 2,4
          1     1     |                        1,1 1,2 1,3 1,4
          0     0     |                                       
        --------------|________________________________________
             ss          -4 -3   -2   -1    0   1   2   3   4 
                       -nx/2     -2   -1    0   1   2      nx/2
           index         0   1    2    3    4   5   6   7   8
                         0                 nx/2            nx"""
        for it in range(self.ntd2,self.nt-1+1):
            for ix in range(1,self.nxd2+1):
                x2[self.nt-it,self.nxd2+ix]=x1[it,ix]
        # Check all coefficients have been overwritten
        x2max=x2.max()
        print('A. x2.max: {0!s}'.format(x2max))
        if x2max>1e19:
            raise UserWarning('Have not successfully overwritten all values.')
        return x2

    def fft2_WK_to_numpy(self,x1):
        """Rearrange 2D FFT from Wheeler-Kiladis format to numpy.fft2 format.

        Inputs:

        <x1> : numpy array, 2D FFT in Wheeler-Kiladis format.

        Output:

        <x2> : numpy array, rearranged version of <x1>, in same format
        that numpy.fft2 delivers.

        """
        # Create empty array for fft2 arrangement, then overwrite
        x2=1e20*np.ones((self.nt,self.nx),dtype=np.cdouble)
        # Convention: it,ix in following code refers to index in *source* array
        # Lower left quadrant. Copy to: 0 to nt/2, 0 to nx/2
        """
        iharm  index  | 
        3,nt/2 3,nt/2 | 3,4 3,3  3,2  3,1  3,0                
          2     2     | 2,4 2,-3 2,-2 2,-1 2,0                
          1     1     | 1,4 1,-3 1,-2 1,-1 1,0                
          0     0     | 0,4 0,3  0,2  0,1  0,0                
        --------------|________________________________________
             ss          -4 -3   -2   -1    0   1   2   3   4 
                       -nx/2     -2   -1    0   1   2      nx/2
           index         0   1    2    3    4   5   6   7   8
                         0                 nx/2            nx

                             copy to

        index
        5,nt-1   |                                        
          4      |                                        
        3,nt/2   | 3,0  3,1  3,2  3,3  3,4               
          2      | 2,0  2,-1 2,-2 2,-3 2,4                 
          1      | 1,0  1,-1 1,-2 1,-3 1,4               
          0      | 0,0  0,1  0,2  0,3  0,4               
                 |________________________________________
                    0    1    2    3    4   5    6     7  x index
                    0                  nx/2           nx-1"""
        for it in range(self.ntd2+1):
            for ix in range(self.nxd2+1):
                x2[it,self.nxd2-ix]=x1[it,ix]
        # Most of upper left quadrant. Copy to: nt/2+1 to nt-1, 1 to nx/2
        """
        iharm  index  | 
        3,nt/2 3,nt/2 | 
          2     2     |                        2,1 2,2 2,3 2,4
          1     1     |                        1,1 1,2 1,3 1,4
          0     0     | 
        --------------|________________________________________
             ss          -4 -3   -2   -1    0   1   2   3   4 
                       -nx/2     -2   -1    0   1   2      nx/2
           index         0   1    2    3    4   5   6   7   8
                         0                 nx/2            nx

                             copy to

        index
        5,nt-1   |      1,1  1,2  1,3  1,4 
          4      |      2,1  2,2  2,3  2,4 
        3,nt/2   |                                       
          2      | 
          1      | 
          0      | 
                 |________________________________________
                    0    1    2    3    4   5    6     7  x index
                    0                  nx/2           nx-1"""
        for it in range(1,self.ntd2-1+1):
            for ix in range(self.nxd2+1,self.nx+1):
                x2[self.nt-it,ix-self.nxd2]=x1[it,ix]
        # Left column of upper left quadrant. Copy to: nt/2+1 to nt-1, 0
        # NB copy the complex conjugate
        """
        iharm  index  | 
        3,nt/2 3,nt/2 | 
          2     2     |                    2,0
          1     1     |                    1,0
          0     0     | 
        --------------|________________________________________
             ss          -4 -3   -2   -1    0   1   2   3   4 
                       -nx/2     -2   -1    0   1   2      nx/2
           index         0   1    2    3    4   5   6   7   8
                         0                 nx/2            nx

                             copy to

        index
        5,nt-1   | 1,0
          4      | 2,0
        3,nt/2   |                                       
          2      | 
          1      | 
          0      | 
                 |________________________________________
                    0    1    2    3    4   5    6     7  x index
                    0                  nx/2           nx-1"""
        for it in range(1,self.ntd2-1+1):
            ix=self.nxd2
            x2[self.nt-it,ix-self.nxd2]=x1[it,ix].conjugate()
        # Lower right quadrant. Copy to: 0 to nt/2, nx/2+1 to nx-1
        # NB copy the complex conjugate
        """
        iharm  index  | 
        3,nt/2 3,nt/2 |                        3,1 3,2 3,3    
          2     2     |                        2,1 2,2 2,3    
          1     1     |                        1,1 1,2 1,3    
          0     0     |                        0,1 0,2 0,3    
        --------------|________________________________________
             ss          -4 -3   -2   -1    0   1   2   3   4 
                       -nx/2     -2   -1    0   1   2      nx/2
           index         0   1    2    3    4   5   6   7   8
                         0                 nx/2            nx

                             copy to

        index
        5,nt-1   | 
          4      | 
        3,nt/2   |                         3,3  3,2   3,1
          2      |                         2,3  2,2   2,1
          1      |                         1,3  1,2   1,1
          0      |                         0,3  0,2   0,1
                 |________________________________________
                    0    1    2    3    4   5    6     7  x index
                    0                  nx/2           nx-1"""
        for it in range(self.ntd2+1):
            for ix in range(self.nxd2+1,self.nx-1+1):
                x2[it,self.nx+self.nxd2-ix]=x1[it,ix].conjugate()
        # Upper right quadrant. Copy to: nt/2+1 to nt-1, nx/2+1 to nx-1
        # NB copy the complex conjugate
        """
        iharm  index  | 
        3,nt/2 3,nt/2 | 
          2     2     |     2,-3 2,-2 2,-1 
          1     1     |     1,-3 1,-2 1,-1 
          0     0     | 
        --------------|________________________________________
             ss          -4 -3   -2   -1    0   1   2   3   4 
                       -nx/2     -2   -1    0   1   2      nx/2
           index         0   1    2    3    4   5   6   7   8
                         0                 nx/2            nx

                             copy to

        index
        5,nt-1   |                         1,-3 1,-2, 1,-1
          4      |                         2,-3 2,-2  2,-1
        3,nt/2   | 
          2      | 
          1      | 
          0      | 
                 |________________________________________
                    0    1    2    3    4   5    6     7  x index
                    0                  nx/2           nx-1"""
        for it in range(1,self.ntd2-1+1):
            for ix in range(1,self.nxd2-1+1):
                x2[self.nt-it,self.nxd2+ix]=x1[it,ix].conjugate()
        x2max=x2.max()
        print('B. x2.max: {0!s}'.format(x2max))
        if x2max>1e19:
            raise UserWarning('Have not successfully overwritten all values.')
        return x2

    def f_filter(self):
        """
        Filter the 2-D Wheeler-Kiladis FFT for selected waves.

        Inputs:

        self.<wave_type> : string describing wave type, from 'none,
        'EK', etc.

        self.<wave_params> : dictionary of various parameters needed to
        selected each wave type, based on wavenumber and
        frequency. Bespoke to each wave type, and designed to mimic
        the dispersion curves of the theoretical waves.

        Outputs:

        Sets attributes self.data_hovfftWKfilt, the filtered 2-D FFT
        (with Wheeler-Kiladis arrangedment of coefficients). Shape
        (nt/2 +1,nx+1).

        """
        # Create array for filter, then overwrite
        x1=np.ones(self.data_hovfftWK.shape,dtype=np.int64)
        # Create wavenumber-frequency filter
        if self.wave_type=='none':
            print('No filtering applied.')
            # Use this to test the whole FFT and reverse FFT process. 
            # As there is no
            # filtering applied, should get back the original Hovmoller.
        elif self.wave_type in ['EK','ER']:
            print('Filter for equatorial waves.')
            # Calculate 2-D fields of freq, ss, omegaf, kk, cphasex for
            # use in creating filters
            freq_coord=self.data_hovfftWK.coord(var_name2long_name['freq'])
            freq_vals=freq_coord.points
            nfreq=nomega=freq_vals.size
            #
            ss_coord=self.data_hovfftWK.coord(var_name2long_name['ss'])
            ss_vals=ss_coord.points
            nss=nkk=ss_vals.size
            print('nfreq,nss: {0!s},{1!s}'.format(nfreq,nss))
            #
            omegaf_coord=self.data_hovfftWK.coord(var_name2long_name['omegaf'])
            omegaf_vals=omegaf_coord.points
            #
            kk_coord=self.data_hovfftWK.coord(var_name2long_name['kk'])
            kk_vals=kk_coord.points
            #
            freq_vals=freq_vals.reshape((nfreq,1))
            ones=np.ones((1,nss))
            self.freq2d=np.dot(freq_vals,ones)
            self.freq2d=create_cube(self.freq2d,self.data_hovfftWK)
            self.freq2d.rename(freq_coord.name())
            self.freq2d.var_name=freq_coord.var_name
            self.freq2d.units=freq_coord.units
            #
            omegaf_vals=omegaf_vals.reshape((nfreq,1))
            ones=np.ones((1,nss))
            self.omegaf2d=np.dot(omegaf_vals,ones)
            self.omegaf2d=create_cube(self.omegaf2d,self.data_hovfftWK)
            self.omegaf2d.rename(omegaf_coord.name())
            self.omegaf2d.var_name=omegaf_coord.var_name
            self.omegaf2d.units=omegaf_coord.units
            #
            ones=np.ones((nfreq,1))
            ss_vals=ss_vals.reshape((1,nss))
            self.ss2d=np.dot(ones,ss_vals)
            self.ss2d=create_cube(self.ss2d,self.data_hovfftWK)
            self.ss2d.rename(ss_coord.name())
            self.ss2d.var_name=ss_coord.var_name
            self.ss2d.units=ss_coord.units
            #
            ones=np.ones((nfreq,1))
            kk_vals=kk_vals.reshape((1,nss))
            self.kk2d=np.dot(ones,kk_vals)
            self.kk2d=create_cube(self.kk2d,self.data_hovfftWK)
            self.kk2d.rename(kk_coord.name())
            self.kk2d.var_name=kk_coord.var_name
            self.kk2d.units=kk_coord.units
            #
            self.cphasex2d=np.where(np.less(np.abs(self.kk2d.data),1e-8),0,self.omegaf2d.data/self.kk2d.data)
            self.cphasex2d=create_cube(self.cphasex2d,self.data_hovfftWK)
            var_name='cphasex'
            self.cphasex2d.rename(var_name2long_name[var_name])
            self.cphasex2d.var_name=var_name
            if self.omegaf2d.units=='rad s-1' and self.kk2d.units=='m-1':
                self.cphasex2d.units='m s-1'
            else:
                raise UserWarning('omegaf or kk have incorrect units.')
            print('cphasex2d. units. max: {0!s},{1!s}'.format(self.cphasex2d.units,self.cphasex2d.data.max()))
            if self.wave_type=='EK':
                print('Filter for equatorial Kelvin waves.')
                print('x1.sum: {0!s}'.format(x1.sum()))
                # Mask out below ss_min
                x1=np.where(np.less(self.ss2d.data,self.wave_params['ss_min']),0,x1)
                print('x1.sum: {0!s}'.format(x1.sum()))
                # Mask out above ss_max
                x1=np.where(np.greater(self.ss2d.data,self.wave_params['ss_max']),0,x1)
                print('x1.sum: {0!s}'.format(x1.sum()))
                # Mask out below freq_min
                if str(self.freq2d.units)!=self.wave_params['freq_units']:
                    raise UserWarning('freq2d and freq parameters must have same units.')
                x1=np.where(np.less(self.freq2d.data,self.wave_params['freq_min']),0,x1)
                print('x1.sum: {0!s}'.format(x1.sum()))
                # Mask out above freq_max
                x1=np.where(np.greater(self.freq2d.data,self.wave_params['freq_max']),0,x1)
                print('x1.sum: {0!s}'.format(x1.sum()))
                # Mask out below cphasex_min
                x1=np.where(np.less(self.cphasex2d.data,self.wave_params['cphasex_min']),0,x1)
                print('x1.sum: {0!s}'.format(x1.sum()))
                # Mask out above cphasex_max
                x1=np.where(np.greater(self.cphasex2d.data,self.wave_params['cphasex_max']),0,x1)
                print('x1.sum: {0!s}'.format(x1.sum()))
            elif self.wave_type=='ER':
                print('Filter for equatorial Rossby waves.')
                raise ToDoError('Code up for ER waves')
            # Create iris cube of filter mask
            self.data_hovfftWKmask=create_cube(x1,self.data_hovfftWK)
        else:
            raise UserWarning('Invalid wave_type.')
        # Apply wavenumber-frequency filter
        self.data_hovfftWKfilt=self.data_hovfftWK*x1

    def f_ifft(self):
        """Calculate inverse 2-D FFT of filtered 2-D Hovmoller.

        Input is the attribute self.data_hovfftWKfilt has the
        Wheeler-Kiladis arrangement of coefficients, and needs to be
        rearranged to the form expected by the numpy.fft routines.

        Then do the inverse 2-D FFT.

        Outputs. Sets attributes:

        self.data_hovfftfilt. This is self.data_hovfftWKfilt with
        coefficients rearranged to the form expected by the numpy.fft
        routines. Shape (nt,nx)

        """
        # Rearrange filtered 2D FFT from Wheeler-Kiladis arrangment to 
        # numpy.fft arrangement
        self.data_hovfftfilt=self.fft2_WK_to_numpy(self.data_hovfftWKfilt.data)
        if self.wave_type=='none' and True:
            print('Checking data_hovfft and data_hovfftfilt are identical.')
            print('Not matching: it,ix,data_hovfft,data_hovfftfilt')
            kount=0
            tol=1
            for it in range(self.nt):
                for ix in range(self.nx):
                    y1=self.data_hovfft[it,ix]
                    y2=self.data_hovfftfilt[it,ix]
                    y3=abs(y1-y2)
                    if y3>tol:
                        print(it,ix,y1,y2)
                        kount+=1
            print('kount: {0!s}'.format(kount))
            if kount>0:
                raise UserWarning('Ouput and input FFTs not identical.')
        # Calculate 2D inverse FFT
        #x2=np.fft.ifft2(self.data_hovfftfilt,norm='ortho')
        x2=np.fft.ifft2(self.data_hovfftfilt)
        # Have checked that the reverse transform of the unfiltered FFT does
        # reproduce the original Hovmoller, to machine accuracy.
        # Note, the output of ifft2 is a complex array, but the imaginary
        # parts are all zero to within machine precision.
        # Hence, take the real part here 
        print('Check that imaginary parts are (near) zero.')
        realmin=x2.real.min()
        realmax=x2.real.max()
        imagmin=x2.imag.min()
        imagmax=x2.imag.max()
        print('realmin,realmax: {0!s}, {1!s}'.format(realmin,realmax))
        print('imagmin,imagmax: {0!s}, {1!s}'.format(imagmin,imagmax))
        tol=1e-8*max(abs(realmin),abs(realmax))
        print('tol: {0!s}'.format(tol))
        if max(abs(imagmin),abs(imagmax))>tol:
            raise UserWarning('Non-zero imaginary parts.')
        x3=x2.real
        #
        # Create iris cube
        x4=create_cube(x3,self.data_hov)
        # Add cell method to describe wavenumber filtering
        cm=iris.coords.CellMethod('point',coords=['time','longitude'],comments='wavenumber filtering. Type: {0.wave_type!s}. Parameters: {0.wave_params!s}'.format(self))
        x4.add_cell_method(cm)
        # Save filtered Hovmoller
        self.data_hovWKfilt=x4
        iris.save(self.data_hovWKfilt,self.file_hovWKfilt)
        if self.archive:
            archive_file(self,self.file_hovWKfilt)
        if self.verbose==2:
            print('write_cube: {0.file_hovWKfilt!s}'.format(self))
        if self.wave_type=='none':
            print('No filtering. Calculating difference between input and output Hovmoller to check code (should be identical to machine precision).')
            self.data_hovdiff=self.data_hovWKfilt-self.data_hov
            xmin=self.data_hov.data.min()
            xmax=self.data_hov.data.max()
            hovscale=max(abs(xmin),abs(xmax))
            tol=1e-6*hovscale
            print('tol: {0!s}'.format(tol))
            print('data_hov.min: {0!s}'.format(xmin))
            print('data_hov.max: {0!s}'.format(xmax))
            xmin=self.data_hovWKfilt.data.min()
            xmax=self.data_hovWKfilt.data.max()
            print('data_hovWKfilt.min: {0!s}'.format(xmin))
            print('data_hovWKfilt.max: {0!s}'.format(xmax))
            xmin=self.data_hovdiff.data.min()
            xmax=self.data_hovdiff.data.max()
            hovdiffscale=max(abs(xmin),abs(xmax))
            print('data_hovdiff.min: {0!s}'.format(xmin))
            print('data_hovdiff.max: {0!s}'.format(xmax))
            if hovdiffscale>tol:
                raise UserWarning('Difference is larger than zero to machine precision.')

#==========================================================================

class CCEWLagrangian(object):

    """Lagrangian database of convectively coupled equatorial waves

    Follows methodology of Baranowski et al. (2016).

    Baranowski DB, Flatau MK, Flatau PJ, Matthews AJ, 2016: Impact of
    atmospheric convectively-coupled Kelvin waves on upper ocean
    variability. J. Geophys. Res., 121, 2045-2059, doi:
    10.1002/2015JD024150.
    
    Called from ccew_lagrangian.py.

    Selected attributes:

    self.wave_type : string denoting type of equatorial wave the
    Lagrangian data base will be calculated for, e.g., 'EK' for
    equatorial Kelvin wave.

    """

    def __init__(self,**descriptor):
        """Initialise from descriptor dictionary.

        Compulsory keywords: 'verbose','source','var_name','level',
        'basedir','filepre', 'time1', 'time2', 'lat1', 'lat2'.
        """
        self.__dict__.update(descriptor)
        self.descriptor=descriptor
        self.name=var_name2long_name[self.var_name]
        source_info(self)
        if self.wave_type in ['EK']:
            self.propagation_direction='eastwards'
        elif self.wave_type in ['ER']:
            self.propagation_direction='westwards'
        else:
            raise UserWarning('Invalid wave_type.')
        # Input files
        ss='_lat_'+str(self.lat1)+'_'+str(self.lat2)+'_'+str(self.time1)[:10]+'_'+str(self.time2)[:10]
        self.file_hov=os.path.join(self.basedir,self.source,'processed',self.var_name+'_'+str(self.level)+self.filepre+'_hov'+ss+'.nc')
        self.file_hovWKfilt=os.path.join(self.basedir,self.source,'processed',self.var_name+'_'+str(self.level)+self.filepre+'_hovWKfilt'+self.wave_type+ss+'.nc')
        # Output files
        self.file_traj_pickle=os.path.join(self.basedir,self.source,'processed',self.var_name+'_'+str(self.level)+self.filepre+'_traj'+self.wave_type+ss+'.pkl')
        self.file_hovtraj=os.path.join(self.basedir,self.source,'processed',self.var_name+'_'+str(self.level)+self.filepre+'_hovtraj'+self.wave_type+ss+'.nc')
        if self.verbose:
            print(self)
        
    def __repr__(self):
        return 'CCEWLagrangian({0.descriptor!r},verbose={0.verbose!r})'.format(self)

    def __str__(self):
        if self.verbose==2:
            ss=h1a+'CCEWLagrangian instance \n'+\
                'var_name: {0.var_name!s} \n'+\
                'level: {0.level!s} \n'+\
                'source: {0.source!s} \n'+\
                'time1: {0.time1!s} \n'+\
                'time2: {0.time2!s} \n'+\
                'wave_type: {0.wave_type!s} \n'+\
                'file_hov: {0.file_hov!s} \n'+\
                'file_hovWKfilt: {0.file_hovWKfilt!s} \n'+\
                'file_traj_pickle: {0.file_traj_pickle!s} \n'+\
                'file_hovtraj: {0.file_hovtraj!s} \n'+h1b
            return ss.format(self)
        else:
            return 'CCEWLagrangian of '+self.source+' '+self.var_name+str(self.level)

    def f_read_hovWKfilt(self):
        """Read WK filtered 2-D time,longitude Hovmoller.

        Creates attributes:

        self.data_hovWKfilt : input 2-D time,longitude Hovmoller of
        equatorial wave-filtered data.

        """
        self.data_hovWKfilt=iris.load_cube(self.file_hovWKfilt)
        # Check it is a 2D (time,longitude) Hovmoller
        self.time_coord=self.data_hovWKfilt.coord('time')
        self.lon_coord=self.data_hovWKfilt.coord('longitude')
        self.ntime=len(self.time_coord.points)
        self.nlon=len(self.lon_coord.points)
        print('ntime,nlon: {0.ntime!s},{0.nlon!s}'.format(self))
        if self.data_hovWKfilt.shape!=(self.ntime,self.nlon):
            raise UserWarning('hovWKfilt does not have correct dimensions.')
        # Also read unfiltered Hovmoller for later use in pruning trajectories
        self.data_hov=iris.load_cube(self.file_hov)
        if self.data_hov.shape!=(self.ntime,self.nlon):
            raise UserWarning('hov does not have correct dimensions.')

    def f_events(self,event_params):
        """Create Lagrangian database of CCEW events.

        Argument:

        event_params : dictionary of parameters for construction of
        Lagrangian data base of CCEW events. These include:

           'threshold' : float, threshold data value (typically
           precipitation rate) which must be exceeded by local maxima
           to be included in data_hovmax

           'threshold_units' : string, e.g., 'mm hr-1', of units of
           data

           'traj_min_time_length' : float, minimum length of time that
           a trajectory must last to be included in the data base,
           e.g., 12.

           'traj_min_length_units' : string, units of
           traj_min_time_length, e.g., 'h'. Must match units of
           self.frequency (i.e., units of time resolution of data
           source).

           'prune_threshold1' and 'prune_threshold2' : floats, used in
           pruning western and eastern ends of raw trajectories to
           take account of over-extension of trajectories based on
           wavenumber-frequency filtered data. More details in
           comments within the function.

        Creates attribute <self.trajectories>, a dictionary of
        dictionaries containing the Lagrangian trajectories. Keys of
        trajectories are integers from 0 upwards.

        Each entry trajectories[index] is a dictionary of information
        on that trajectory. Keys are:

           'npts' : an integer, the number of (time,longitude) points
           on the trajectory.

           'time_indices' and 'lon_indices' : these are two lists (of
           the length npts and to be paired together) of the time and
           longitude indices of the maxima (the 1's) in
           self.data_hovmax.

           'times' : a list of length npts of datetime.datetime or
           similar objects corresponding to the actual times of the
           time_indices.

           'lons' : a list of length npts of longitude values
           corresponding to the actual longitudes of the lon_indices.

           'tunits' : a CF units object of the times units.

           'lonunits': a CF units object of the longitude units.

           'data_vals' : a list of length npts of the data values at
           the points along the trajectory from the WK-filtered
           Hovmoller self.data_hovWKfilt.

           'dataunits' : a CF units object of the data units (e.g., mm
           hr-1).

        """
        def clamp(n,minn,maxn):
            """Clamp n within range minn<=n<=maxn."""
            if n<minn:
                return minn
            elif n>maxn:
                return maxn
            else:
                return n
        self.event_params=event_params
        # Create self.data_hovmax
        # Check that units of threshold are same as units of input Hovmoller
        if self.event_params['threshold_units']!=str(self.data_hovWKfilt.units):
            raise UserWarning('threshold and Hovmoller must have same units.')
        # First, create an array of zeroes same size as input Hovmoller
        x1=np.zeros(self.data_hovWKfilt.shape,dtype=np.int64)
        # At each longitude, find local maxima in time (by comparison with
        # neighbours) and overwrite with a '1'
        for ilonc in range(self.nlon):
            x2=self.data_hovWKfilt.data[:,ilonc]
            # Find indices of all maxima
            x3=scipy.signal.argrelmax(x2)
            #print('x3: {0!s}'.format(x3))
            # Then check that data value is above threshold at these maxima
            # Note scipy.signal.find_peaks is not available in scipy!
            for time_index in x3[0]:
                x2a=x2[time_index]
                if x2a>=self.event_params['threshold']:
                    x1[time_index,ilonc]=1
                    #print('time_index,x2a,True: {0!s},{1!s}'.format(time_index,x2a))
                else:
                    #print('time_index,x2a,False: {0!s},{1!s}'.format(time_index,x2a))
                    pass
        self.data_hovmax=create_cube(x1,self.data_hovWKfilt)
        var_name='dummy'
        long_name=var_name2long_name[var_name]
        self.data_hovmax.rename(long_name)
        self.data_hovmax.var_name=var_name
        self.data_hovmax.units='1'
        # Create copy of data_hovmax to calculate Lagrangian trajectories from
        # hovmax1: once a trajectory has been identified, overwrite its 1s
        #  with 0s, and search it again for the next trajectory
        hovmax1=copy.deepcopy(self.data_hovmax)
        xx=hovmax1.data.sum()
        print('hovmax1.data.sum: {0!s}'.format(xx))
        tcoord=self.data_hovmax.coord('time')
        tunits=tcoord.units
        loncoord=self.data_hovmax.coord('longitude')
        lonunits=loncoord.units
        dataunits=self.data_hovWKfilt.units
        # Initialise
        self.trajectories={}
        kevent=0
        ktraj_include=0
        # Set minimum allowed lifetime of trajectory to be included in
        # data base (directory of directories)
        # Although this was a requirement in Baranowski et al. (2016a)
        # and is implemented here, in practice it never (or almost
        # never) excludes any trajectories.
        if self.frequency[-1]==self.event_params['traj_min_time_length_units']:
            # Minimum allowed lifetime of trajectory is 
            # self.event_params['traj_min_time_length']
            # Convert this to number of time indices that a trajectory
            # must cover to be accepted
            self.traj_min_time_length_indices=divmod(self.event_params['traj_min_time_length'],int(self.frequency[:-1]))[0]+1
            print('traj_min_time_length_indices: {0!s}'.format(self.traj_min_time_length_indices))
        else:
            raise UserWarning('Units of frequency of input data and min lifetime must match.')
        # Set nevent_max to more events than you expect to find. nevent_max
        # is just there to stop an endless loop in case detection of 
        # reaching last time in data fails.
        nevent_max=10000
        if self.propagation_direction=='eastwards':
            # Search through copy of data_hovmax to find trajectories.
            # Start at 'bottom left', i.e., first time and furthest westward
            # point, then work eastwards in longitude, then forward in time.
            # Once the beginning of a trajectory has been found, transcribe
            # its propagation path to the Lagrangian data base, and overwrite
            # it in the (copy of the) Eulerian data base with zeros. Then 
            # repeat until all the trajectories have been found.
            reached_last_time=False
            time_index_start=0
            while kevent<nevent_max and not reached_last_time:
                print('### CCEW kevent,time_index_start: {0!s},{1!s}'.format(kevent,time_index_start))
                for time_index in range(time_index_start,self.ntime):
                    for lon_index in range(self.nlon):
                        x2=hovmax1.data[time_index,lon_index]
                        #print('time_index,lon_index,x2: {0!s},{1!s},{2!s}'.format(time_index,lon_index,x2))
                        if x2==1:
                            # Found a 1
                            break
                    if x2==1:
                        break
                print('Start of new trajectory time_index,lon_index {0!s},{1!s}'.format(time_index,lon_index))
                # Build up entry for this trajectory
                time_indices=[time_index,]
                lon_indices=[lon_index,]
                finding_current_traj_points=True
                while finding_current_traj_points:
                    # Check that following time point at current longitude is 0
                    time_indexp1=time_index+1
                    if time_indexp1>=self.ntime:
                        print('Reached final time. Stopping.')
                        reached_last_time=True
                        break
                    if hovmax1.data[time_indexp1,lon_index]!=0:
                        raise UserWarning('Two consecutive in time 1s found at same longitude.')
                    # Find the following 1 in the grid point to the east of
                    # the current grid point.
                    # It must be found in the following 6 hours
                    lon_index+=1
                    if lon_index==self.nlon:
                        lon_index=0 # longitude is periodic
                    if self.frequency=='3h':
                        n_time_check=3 # 3 timesteps covers 6 hours (t=0,3,6 hr)
                    else:
                        raise ToDoError('Code up for other time resolutions.')
                    found_a_1=False
                    #print('Start of loop at time_index,lon_index: {0!s},{1!s}'.format(time_index,lon_index))
                    for i_time_check in range(n_time_check):
                        time_index_check=time_index+i_time_check
                        if time_index_check>=self.ntime:
                            print('Reached final time. Stopping.')
                            reached_last_time=True
                            break
                        x2=hovmax1.data[time_index_check,lon_index]
                        #print('i_time_check,time_index_check,lon_index,x2: {0!s},{1!s},{2!s},{3!s}'.format(i_time_check,time_index_check,lon_index,x2))
                        if x2==1:
                            found_a_1=True
                            time_indices.append(time_index_check)
                            lon_indices.append(lon_index)
                            time_index=time_index_check
                            #print('Trajectory point at i_time_check,time_index,lon_index: {0!s},{1!s},{2!s}.'.format(i_time_check,time_index,lon_index))
                            break
                    if not found_a_1:
                        # Reached end of i_time_check loop and not found a 1.
                        # Trajectory ends
                        finding_current_traj_points=False
                #print('Final trajectory: time_indices,lon_indices: {0!s},{1!s}'.format(time_indices,lon_indices))
                # If trajectory lifetime is greater than minimum allowed
                # add to data base
                trajc_time_length_indices=time_indices[-1]-time_indices[0]+1
                if not reached_last_time:
                    if trajc_time_length_indices>=self.traj_min_time_length_indices:
                        print('Accepting trajectory ktraj_include: {0!s}.'.format(ktraj_include))
                        trajc={}
                        trajc['time_indices']=time_indices
                        trajc['lon_indices']=lon_indices
                        self.trajectories[ktraj_include]=trajc
                        ktraj_include+=1
                    else:
                        print('######## Rejecting current trajectory, too short.')
                # Overwrite the 1s from the current trajectory with 0s in
                # hovmax1, ready to search hovmax1 for the next trajectory
                # NB hovmax1.data.sum() is the number of 1s (left) in hovmax
                # This will decrease as each trajectory is removed and
                # should be zero at the end
                nind=len(time_indices)
                for ii in range(nind):
                    hovmax1.data[time_indices[ii],lon_indices[ii]]=0
                xx=hovmax1.data.sum()
                print('nind,hovmax1.data.sum: {0!s},{1!s}'.format(nind,xx))
                kevent+=1
                # Reset time_index_start so do not have to go through
                # entire hovmax for each trajectory
                time_index_start=time_indices[0] 
            if kevent>=nevent_max:
                raise UserWarning('Need to increase nevent_max.')
            if xx!=0:
                ss="""####################################################

                There should be no 1s left in hovmax at the end.

                ######################################################"""
                print(ss)
            if xx>300: # Arbitrary cutoff
                raise UserWarning('Sorry, that is too many 1s left.')
        elif self.propagation_direction=='westwards':
            raise ToDoError('Code up starting from "right hand" end and working leftwards! Best to extend code above to do this rather than repeat.')
        # Prune western and eastern ends of trajectories using unfiltered
        # data to avoid over extension of trajectories because of wavenumber-
        # frequency filtering
        for keyc in self.trajectories.keys():
            print('Pruning keyc: {0!s}'.format(keyc))
            time_indices=self.trajectories[keyc]['time_indices']
            npts=len(time_indices)
            lon_indices=self.trajectories[keyc]['lon_indices']
            data_unfiltered_vals=[self.data_hov.data[time_indices[xx],lon_indices[xx]] for xx in range(npts)]
            data_vals=[self.data_hovWKfilt.data[time_indices[xx],lon_indices[xx]] for xx in range(npts)]
            # Start at western end of trajectory and check that 9 hr
            # mean exceeds threshold of prune_threshold1 (e.g., 0.5 mm hr-1)
            # Then check that this 9 hr mean, averaged over all points
            # on the trajectory from the current point to the
            # (eastward) point 12 hr later exceeds a threshold of
            # prune_threshold2 (e.g., 0.25 mm hr-1).
            # If both criteria are not met, prune this point from the
            # (westward) end of the trajectory.
            # Move to the next (eastward) point and repeat, until
            # the criteria is satisfied.
            if self.frequency=='3h':
                data_unfiltered_plus3h_vals=[self.data_hov.data[clamp(time_indices[xx]+1,0,self.ntime),lon_indices[xx]] for xx in range(npts)]
                data_unfiltered_minus3h_vals=[self.data_hov.data[clamp(time_indices[xx]-1,0,self.ntime),lon_indices[xx]] for xx in range(npts)]
                data_unfiltered_9h_mean=[(data_unfiltered_minus3h_vals[xx]+data_unfiltered_vals[xx]+data_unfiltered_plus3h_vals[xx])/3.0 for xx in range(npts)]
                k12h=5 # for 3 hr data
            else:
                raise ToDoError('Code up for non-3 hourly data.')
            threshold1=self.event_params['prune_threshold1']
            threshold2=self.event_params['prune_threshold2']
            print('threshold1,threshold2: {0!s},{1!s}'.format(threshold1,threshold2))
            # First, find the indices of the points to be pruned
            prune_western_indices=[]
            for xx in range(npts):
                if data_unfiltered_9h_mean[xx]>=threshold1:
                    criterium1=True
                else:
                    criterium1=False
                if sum(data_unfiltered_9h_mean[xx:xx+k12h])/float(k12h)>=threshold2:
                    criterium2=True
                else:
                    criterium2=False
                #print('xx,criterium1,criterium2: {0!s},{1!s},{2!s}'.format(xx,criterium1,criterium2))
                if not(criterium1 and criterium2):
                    # Both criteria not satisfied. Mark for pruning.
                    prune_western_indices.append(xx)
                else:
                    # Both criteria satisfied. Stop process.
                    break
            print('prune_western_indices: {0!s}'.format(prune_western_indices))
            # Repeat from eastern end of trajectory, moving westwards
            prune_eastern_indices=[]
            for xx in range(npts-1,0-1,-1):
                if data_unfiltered_9h_mean[xx]>=threshold1:
                    criterium1=True
                else:
                    criterium1=False
                if sum(data_unfiltered_9h_mean[xx-k12h+1:xx+1])/float(k12h)>=threshold2:
                    criterium2=True
                else:
                    criterium2=False
                #print('xx,criterium1,criterium2: {0!s},{1!s},{2!s}'.format(xx,criterium1,criterium2))
                if not(criterium1 and criterium2):
                    # Both criteria not satisfied. Mark for pruning.
                    prune_eastern_indices.append(xx)
                else:
                    # Both criteria satisfied. Stop process.
                    break
            print('prune_eastern_indices: {0!s}'.format(prune_eastern_indices))
            # Now do the actual pruning
            # Do this by retaining the non-pruned points
            if len(prune_western_indices)>0:
                x1=prune_western_indices[-1]
            else:
                x1=-1
            if len(prune_eastern_indices)>0:
                x2=prune_eastern_indices[-1]
            else:
                x2=npts
            print('x1,x2,npts: {0!s},{1!s},{2!s}'.format(x1,x2,npts))
            self.trajectories[keyc]['time_indices']=time_indices[x1+1:x2]
            self.trajectories[keyc]['lon_indices']=lon_indices[x1+1:x2]
        # Now rerun check on trajectory lifetime after pruning
        # This check is likely to remove a significant fraction of the
        # trajectories.
        xx1={}
        ktraj_include=0
        for keyc in self.trajectories.keys():
            time_indices=self.trajectories[keyc]['time_indices']
            npts=len(time_indices)
            lon_indices=self.trajectories[keyc]['lon_indices']
            if npts==0:
                trajc_time_length_indices=0
            else:
                trajc_time_length_indices=time_indices[-1]-time_indices[0]+1
            print('Checking trajectory length again: keyc,npts: {0!s},{1!s}'.format(keyc,npts))
            if trajc_time_length_indices>=self.traj_min_time_length_indices:
                print('Accepting trajectory ktraj_include: {0!s}.'.format(ktraj_include))
                trajc={}
                trajc['time_indices']=time_indices
                trajc['lon_indices']=lon_indices
                xx1[ktraj_include]=trajc
                ktraj_include+=1
            else:
                print('######## Rejecting current trajectory, too short.')
        self.trajectories=xx1
        # Add time stamps and actual longitudes (in addition to time 
        # and longitude indices) to trajectories dictionary
        # Also add data values along trajectories
        for keyc in self.trajectories.keys():
            time_indices=self.trajectories[keyc]['time_indices']
            npts=len(time_indices)
            times=[tunits.num2date(tcoord.points[xx]) for xx in time_indices]
            lon_indices=self.trajectories[keyc]['lon_indices']
            lons=[loncoord.points[xx] for xx in lon_indices]
            data_vals=[self.data_hovWKfilt.data[time_indices[xx],lon_indices[xx]] for xx in range(npts)]
            self.trajectories[keyc]['npts']=npts
            self.trajectories[keyc]['times']=times
            self.trajectories[keyc]['lons']=lons
            self.trajectories[keyc]['lonunits']=lonunits
            self.trajectories[keyc]['data_vals']=data_vals
            self.trajectories[keyc]['dataunits']=dataunits
        # Save trajectories as pickle file
        # NB to re-read from the pickle file, use:
        # trajectories=pickle.load(filec,'rb')
        filec=open(self.file_traj_pickle,'wb')
        pickle.dump(self.trajectories,filec)
        filec.close()
        # Create another hovmoller of 0s and 1s to represent the final
        # trajectories. This will be a modification of self.data_hovmax
        # to take acount of points lost through pruning and rejection of
        # short trajectories
        xx1=np.zeros(self.data_hovmax.data.shape,dtype=np.int64)
        for keyc in self.trajectories.keys():
            time_indices=self.trajectories[keyc]['time_indices']
            npts=len(time_indices)
            lon_indices=self.trajectories[keyc]['lon_indices']
            for ii in range(npts):
                xx1[time_indices[ii],lon_indices[ii]]=1
        self.data_hovtraj=create_cube(xx1,self.data_hovmax)
        var_name='dummy'
        long_name=var_name2long_name[var_name]
        self.data_hovtraj.rename(long_name)
        self.data_hovtraj.var_name=var_name
        self.data_hovtraj.units='1'
        xx=self.data_hovmax.data.sum()
        print('data_hovmax.data.sum: {0!s}'.format(xx))
        xx=self.data_hovtraj.data.sum()
        print('data_hovtraj.data.sum: {0!s}'.format(xx))
        iris.save(self.data_hovtraj,self.file_hovtraj)
        if self.archive:
            archive_file(self,self.file_hovtraj)

    def f_read_lagrangian(self):
        """Read CCEW Lagrangian data base

        Creates attributes:

        self.trajectories : as defined in f_events()

        """
        filec=open(self.file_traj_pickle,'rb')
        self.trajectories=pickle.load(filec)
        filec.close()
        if self.verbose==2:
            for keyc in self.trajectories:
                print('keyc,start_time,start_lon: {0!s},{1!s},{2!s}'.format(keyc,self.trajectories[keyc]['times'][0],self.trajectories[keyc]['lons'][0]))

    def f_create_time_domain(self,tdomain_params):
        """Create time domain from Lagrangian data base of CCEWs.

        Arguments:

        tdomain_params : Dictionary of parameters used to select dates
        from the CCEW Lagrangian event data base, that then go into
        the time domain. Keys are:

            'lonc' : float or int, longitude that the CCEW events must
            cross, e.g., 90.0. The time that goes into the time domain
            for each event is the time the CCEW crosses longitude
            lonc. Value of lonc must be set.

            'threshold' : float or False, optional threshold value
            that must be exceeded by CCEW as it passes through lonc,
            e.g., 0.2. Set to False to disable.

            'threshold_units' : string, units of data for threshold,
            e.g., 'mm hr-1'.

            'min_lon_extent' : float or False, optional minimum
            longitudinal extent that the CCEW event must extend over,
            e.g., 20.0. This is split equally either side of lonc. For
            example, if lonc is 90, and min_lon_extent is 20, then to
            be included in the time domain, the CCEW event must extend
            between at least 80 and 100 degrees east. Set to False to
            disable.

            'round_to_nearest_time' : string of False. If False, pass
            datetime of crossing lonc to time domain at full
            resolution, e.g., 3 hour resolution. If 'd', then round
            datetime to nearest day at 00 UTC. If '6h', then round
            datetime to nearest 00, 06, 12, 18 UTC.

        """
        self.tdomain_params=tdomain_params
        self.lonc=self.tdomain_params['lonc']
        self.threshold=self.tdomain_params['threshold']
        self.threshold_units=self.tdomain_params['threshold_units']
        self.min_lon_extent=self.tdomain_params['min_lon_extent']
        self.round_to_nearest_time=self.tdomain_params['round_to_nearest_time']
        # Create time domain idx (name) and headers
        # Permission to be ad hoc here!
        idx='CC'+self.wave_type+str(self.lonc)+'E'+str(self.time1.year)[-2:]+'-'+str(self.time2.year)[-2:]
        if self.threshold:
            idx=idx+'-'+str(self.threshold)
        if self.round_to_nearest_time=='d':
            idx=idx+'-00UTC'
        elif self.round_to_nearest_time=='6h':
            idx=idx+'-6h'
        elif not(self.round_to_nearest_time):
            idx=idx+'-'+self.frequency
        else:
            raise UserWarning('Invalid round_to_nearest_time.')
        print('idx: {0!s}'.format(idx))
        header1='# CC'+self.wave_type+'W arrival times at basepoint.'
        if self.round_to_nearest_time:
            header1=header1+' Rounded to nearest '+self.round_to_nearest_time
        header1=header1+'\n'
        header2='# Basepoint '+str(self.lonc)+'E. '+str(self.time1.year)+'-'+str(self.time2.year)+'.'
        if self.threshold:
            header2=header2+' Threshold '+str(self.threshold)+' '+self.threshold_units+'.'
        header2=header2+'\n'
        print('header1: {0!s}'.format(header1))
        print('header2: {0!s}'.format(header2))
        # Extract crossing times from Lagrangian data base
        # (Apply extra thresholds if necessary here)
        crossing_times=[]
        for keyc in self.trajectories.keys():
            lons=self.trajectories[keyc]['lons']
            lon1=lons[0]
            lon2=lons[-1]
            times=self.trajectories[keyc]['times']
            time1=times[0]
            time2=times[-1]
            npts=self.trajectories[keyc]['npts']
            data_vals=self.trajectories[keyc]['data_vals']
            if lon1<=self.lonc<=lon2:
                timec=time1
                for ii in range(npts):
                    if lons[ii]>=self.lonc:
                        break
                    # ii is now the index of the crossing point
                    timec=times[ii]
                # Apply amplitude threshold at crossing longitude if required
                amp_flag=True
                if self.threshold:
                    if self.threshold_units!=self.trajectories[keyc]['dataunits']:
                        raise UserWarning('Units of threshold and data must match.')
                    if data_vals[ii]<self.threshold:
                        amp_flag=False
                # Round time if required
                if self.round_to_nearest_time=='d':
                    if self.calendar=='gregorian':
                        xx=datetime.datetime(timec.year,timec.month,timec.day)
                        if timec.hour>=12:
                            xx=xx+datetime.timedelta(days=1)
                        timec=xx
                    else:
                        raise ToDoError('Code up for non-Gregorian calendar.')
                elif self.round_to_nearest_time=='6h':
                    if self.calendar=='gregorian':
                        deltahour=6
                        x1,x2=divmod(timec.hour,deltahour)
                        xx=datetime.datetime(timec.year,timec.month,timec.day,x1*deltahour) # round down to previous 6h
                        if x2>=deltahour/2:
                            xx=xx+datetime.timedelta(hours=deltahour) # round up to next 6h
                        print('timec,xx: {0!s}, {1!s}'.format(timec,xx))
                        timec=xx
                    else:
                        raise ToDoError('Code up for non-Gregorian calendar.')
                if amp_flag:
                    print('Including: keyc,lon1,lonc,lon2,time1,timec,time2: {0!s}, {1!s}, {2!s}, {3!s}, {4!s}, {5!s}, {6!s}'.format(keyc,lon1,self.lonc,lon2,time1,timec,time2))
                    crossing_times.append(timec)
        # Sort by time as this is not guaranteed by the methodology
        crossing_times.sort()
        # Create list of paired (start_time,end_time) at (00UTC,23:59:59UTC)
        if self.calendar=='gregorian':
            crossing_times2=[[datetime.datetime(xx.year,xx.month,xx.day),datetime.datetime(xx.year,xx.month,xx.day,23,59,59)] for xx in crossing_times]
        else:
            raise ToDoError('Code up for non-Gregorian calendar.')
        # Create time domain
        tdomain=TimeDomain(idx)
        tdomain.header=(header1,header2)
        tdomain.datetimes=crossing_times # single values
        #tdomain.datetimes=crossing_times2 # paired values
        tdomain.datetime2ascii()
        tdomain.write_ascii()
