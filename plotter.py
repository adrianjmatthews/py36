"""Plot labels, simple functions for use in plotting. """

import numpy as np

import iris
import matplotlib as mpl


months_January={1:'January', 2:'February', 3:'March', 4:'April', 5:'May', 6:'June', 7:'July', 8:'August', 9:'September', 10:'October', 11:'November', 12:'December'}
months_Jan={1:'Jan', 2:'Feb', 3:'Mar', 4:'Apr', 5:'May', 6:'Jun', 7:'Jul', 8:'Aug', 9:'Sep', 10:'Oct', 11:'Nov', 12:'Dec'}
months_1_Jan={1:'1 Jan', 2:'1 Feb', 3:'1 Mar', 4:'1 Apr', 5:'1 May', 6:'1 Jun', 7:'1 Jul', 8:'1 Aug', 9:'1 Sep', 10:'1 Oct', 11:'1 Nov', 12:'1 Dec'}
months_J={1:'J', 2:'F', 3:'M', 4:'A', 5:'M', 6:'J', 7:'J', 8:'A', 9:'S', 10:'O', 11:'N', 12:'D'}

def lonlat2string(val,lonlat,ndigits=2,format=0):
    """Return a string of form e.g., '5$^\{circ}$N'.

    Use to create longitude and latitude axes labels for plotting.

    Inputs:

    val is a float or integer, e.g., 5, -112.5

    lonlat is a string, either 'lon' or 'lat'

    format : 0, format strings with degree sign
             1, format strings without degree sign

    Outputs:

    lonlatstring is a string, e.g., '5$^\{circ}$S'.

    """
    # Check lonlat is valid
    if lonlat not in ['lon','lat']:
        raise UserWarning("lonlat not valid. Must be 'lon' or 'lat'.")
    # Remove decimal point if integer value, eg 5.0 becomes 5
    #   and take absolute value
    # If not an integer value, round to ndigits decimal points
    if int(val)==val:
        xx=abs(int(val))
    else:
        xx=abs(round(val,ndigits))
    # Set sign
    if lonlat=='lon':
        if val==-180:
            xx=180
            sign='E'
        if -180<val<0:
            sign='W'
        elif 0<=val<=180:
            sign='E'
        elif 180<val<360:
            xx=360-xx
            sign='W'
        elif val==360:
            xx=0
            sign='E'
        else:
            raise ValueError('longitude must be in range -180 <= val <= 360')
    else:
        if -90<=val<0:
            sign='S'
        elif 0<=val<=90:
            sign='N'
        else:
            raise ValueError('latitude must be in range -90 <= val <= 90')
    # Create string
    if format==0:
        str1='$^\circ$'
    elif format==1:
        str1=''
    else:
        raise UserWarning('Invalid format argument.')
    lonlatstring=str(xx)+str1+sign
    return lonlatstring

def suppress_label(labels,tickmarks,val):
    """Set specific label to emptry string."""
    print('# suppress_label.')
    print('Input labels: {0!s}'.format(labels))
    index=list(tickmarks).index(val)
    labels[index]=''
    print('Output labels: {0!s}'.format(labels))
    return labels

def reverse_colormap(cmap,name='my_cmap_r'):
    """
    In: 
    cmap, name 
    Out:
    my_cmap_r

    Explanation:
    t[0] goes from 0 to 1
    row i:   x  y0  y1 -> t[0] t[1] t[2]
                   /
                  /
    row i+1: x  y0  y1 -> t[n] t[1] t[2]

    so the inverse should do the same:
    row i+1: x  y1  y0 -> 1-t[0] t[2] t[1]
                   /
                  /
    row i:   x  y1  y0 -> 1-t[n] t[2] t[1]
    """        
    reverse=[]
    k=[]   
    for key in cmap._segmentdata:    
        k.append(key)
        channel=cmap._segmentdata[key]
        data=[]
        for t in channel:                    
            data.append((1-t[0],t[2],t[1]))            
        reverse.append(sorted(data))    
    LinearL=dict(zip(k,reverse))
    my_cmap_r=mpl.colors.LinearSegmentedColormap(name,LinearL) 
    return my_cmap_r

def fixed_aspect_ratio(axc,ratio=1):
    """Set the physical aspect ratio of an axis.

    axc is a matplotlib axis, with x and y limits already set.
    
    ratio is the desired x/y ratio.

    Redundant.  Use panel-plots to set location of axes on figure
    """
    #raise DeprecationWarning('Do not use.  Use panel-plots instead.')
    xvals,yvals=axc.get_xlim(),axc.get_ylim()
    xrange=abs(xvals[1]-xvals[0])
    yrange=abs(yvals[1]-yvals[0])
    axc.set_aspect(ratio*(xrange/yrange),adjustable='box')

def cube_reset_time_coord(cube,verbose=True):
    """Reset time coord on iris cube so time begins at zero.

    This helps matplotlib, which does not seem to be able to deal with
    very large time values, eg 1887432 (hours since 1800-01-01).

    Just call this function from plotting scripts, rather than general
    analysis scripts.

    Input is <cube>

    No returned output.  The input <cube> has its time axis changed. 
    """
    # Extract existing time coordinate
    time_coord=cube.coord('time')
    time_units=time_coord.units
    dates=time_units.num2date(time_coord.points)
    time_interval=time_units.name.split()[0] # e.g., 'hour'
    time0=time_units.num2date(time_coord.points[0]) # first time value
    # Create new time coordinate
    new_time_units_str=time_interval+' since '+str(time0)
    new_time_coord_data=time_units.convert(time_coord.points,new_time_units_str)
    new_time_coord=iris.coords.DimCoord(new_time_coord_data,standard_name='time',units=new_time_units_str)
    # Replace time coordinate on cube
    coords=cube.coords()
    icoord=False
    for icoord in range(len(coords)):
        if coords[icoord].name()=='time':
            time_coord_index=icoord
    cube.remove_coord('time')
    cube.add_dim_coord(new_time_coord,time_coord_index)
    if verbose:
        print('Resetting time coordinate: {0!s}'.format(new_time_coord.points))

def conv_time_units(cube1,cube2):
    """Convert time axis of cube1 to have same time units as cube 2.

    This is useful when plotting both cube1 and cube2 on the same
    subplot in matplotlib.
    """
    time_coord1=cube1.coord('time')
    time_units1=time_coord1.units
    #
    time_coord2=cube2.coord('time')
    time_units2=time_coord2.units
    #
    new_time_vals=[time_units2.date2num(time_units1.num2date(xx)) for xx in time_coord1.points]
    new_time_coord=iris.coords.DimCoord(new_time_vals,standard_name='time',units=time_units2)
    #
    coord_names=[dimc.standard_name for dimc in cube1.dim_coords]
    time_index=coord_names.index('time')
    cube1.remove_coord('time')
    cube1.add_dim_coord(new_time_coord,time_index)

def levels_list(start,end,interval,zero=True,ndigits=-1):
    """Return a list of (contour) levels.

    Runs from <start> to <end> with <step> of interval.

    If <zero> is False, removes 0 from the list (if it exists),
    otherwise leaves it.

    If <ndigits> is True and an integer value, all levels are rounded
    to a decimal with ndigits decimal places using the round function.
    Use ndigits=0 for integer values.

    Returns <levels>."""

    print('start={0!s}; end={1!s}; interval={2!s}; zero={3!s}; ndigits={4!s}'.format(start,end,interval,zero,ndigits))
    delta=1e-12
    levels=[]
    for xx in np.arange(start,end+delta,interval):
        levelc=xx
        if ndigits>=0:
            levelc=round(xx,ndigits)
            if ndigits==0:
                levelc=int(levelc)
        if not(levelc==0 and not zero):
            levels.append(levelc)
        #print(xx,levelc)
    print('levels : {0!s}'.format(levels))
    nlevels=len(levels)
    print('nlevels : {0!s}'.format(nlevels))
    return levels

def contour_prune(cube,levels_in,colors_in,linestyles_in,linewidths_in,verbose=True):
    """Prunes contour parameters dependent on range of data to plot.

    The matplotlib.contour method, incredibly, does *not* match up the
    lists of levels, colors, linestyles and linewidths that the user
    supplies. If the minimum value of the data to be plotted is
    greater than the first (lowest) level, this level is not contoured
    (naturally). However, the first value of the colors, linestyles
    and linewidths lists are *not* also discarded. They are reserved,
    to be used with the first level that is actually contoured. 

    In the following example, the user has supplied a list of levels
    and wishes for negative contours to be blue and solid, the zero
    contour to be black, solid with double linewidth, and positive
    contours to be red and dashed.

    levels=[-2,-1,0,1,2]
    colors=['blue','blue','black','red','red']
    linestyles=['-','-','-','--','--']
    linewidths=[1,1,2,1,1]

    If the minimum values of the data to be plotted is < -2, this will
    happen. However, if the minimum value is, e.g., -1.8, the
    following contours will be plotted:

    -1, blue, '-', 1
    0,  blue, '-', 1
    1, black, '-', 2
    2, red,   '--', 1

    This is because matplotlib, incredibly, only starts using the
    supplied colorors etc, once it plots the first level, which may or
    may not be the same level supplied by the user.

    Solution. Call this function immediately before plotting, to prune
    the colors etc., so the contoured levels do match up with the
    desired colors etc.

    Inputs:

    <cube> : 2-D iris cube to be contoured
    <levels_in> : list of levels to be contoured
    <colors_in> : list of colors to be used, with 1:1 correspondence with levels_in, or single value
    <linestyles_in> : ditto for linestyles
    <linewidths_in> : ditto for linewidths

    Outputs:

    <levels_out> : pruned list of levels to be contoured
    <colors_out> : list of colors to be used, with 1:1 correspondence with levels_out
    <linestyles_out> : ditto for linestyles
    <linewidths_out> : ditto for linewidths

    Return levels_out,colors_out,linestyles_out,linewidths_out
    """

    # Check that list of input levels etc. all have same length
    nn=len(levels_in)
    check=True
    if (type(colors_in) is list or type(colors_in) is tuple) and len(colors_in)!=nn: check=False
    if (type(linestyles_in) is list or type(linestyles_in) is tuple) and len(linestyles_in)!=nn: check=False
    if (type(linewidths_in) is list or type(linewidths_in) is tuple) and len(linewidths_in)!=nn: check=False
    if not check:
        raise ValueError('Input list of levels etc. must all have same length.')
    # Check that cube is 2-D
    if len(cube.data.shape)!=2:
        raise ValueError('Cube must be 2-D.')
    # Find minimum and maximum data values to be plotted.
    data_min=cube.data.min()
    data_max=cube.data.max()
    # Prune levels etc. as needed to lie within data_min,data_max range
    levels_out=[]
    colors_out=[]
    linestyles_out=[]
    linewidths_out=[]
    for ii in range(len(levels_in)):
        if data_min<=levels_in[ii]<=data_max:
            levels_out.append(levels_in[ii])
            if type(colors_in) is list or type(colors_in) is tuple: colors_out.append(colors_in[ii])
            if type(linestyles_in) is list or type(linestyles_in) is tuple: linestyles_out.append(linestyles_in[ii])
            if type(linewidths_in) is list or type(linewidths_in) is tuple: linewidths_out.append(linewidths_in[ii])
    # If colors_in are just a single value (not a list), return this input value
    if type(colors_in) is not list and type(colors_in) is not tuple: colors_out=colors_in
    if type(linestyles_in) is not list and type(linestyles_in) is not tuple: linestyles_out=linestyles_in
    if type(linewidths_in) is not list and type(linewidths_in) is not tuple: linewidths_out=linewidths_in
    if verbose:
        print('# contour_prune')
        print('data_min,data_max: {0!s},{1!s}'.format(data_min,data_max))
        print('levels_in: {0!s}'.format(levels_in))
        print('levels_out: {0!s}'.format(levels_out))
        print('colors_in: {0!s}'.format(colors_in))
        print('colors_out: {0!s}'.format(colors_out))
        print('linestyles_in: {0!s}'.format(linestyles_in))
        print('linestyles_out: {0!s}'.format(linestyles_out))
        print('linewidths_in: {0!s}'.format(linewidths_in))
        print('linewidths_out: {0!s}'.format(linewidths_out))
        print('# end contour_prune')
    # Check that list of output levels etc. all have same length
    nn=len(levels_out)
    check=True
    if (type(colors_out) is list or type(colors_out) is tuple) and len(colors_out)!=nn: check=False
    if (type(linestyles_out) is list or type(linestyles_out) is tuple) and len(linestyles_out)!=nn: check=False
    if (type(linewidths_out) is list or type(linewidths_out) is tuple) and len(linewidths_out)!=nn: check=False
    if not check:
        raise ValueError('Output list of levels etc. must all have same length.')

    return levels_out,colors_out,linestyles_out,linewidths_out

def create_regular_polygon(x_centre=0,y_centre=0,radius=1,npts=24,verbose=False):
    """Calculate the x,y coordinates of the vertices of a regular polygon.

    Inputs:

    <x_centre>, <y_centre> : the x,y coordinates of the centre of the polygon.

    <radius> : the radius of a bounding circle that passes through all
    the vertices of the polygon.

    <npts> : the number of vertices of the polygon, e.g., 4 for a square.

    Outputs:

    <xvals>, <yvals> : lists of length <npts> of the x and y coordinates of
    the vertices of the polygon.

    Usage: if <npts> is large, the polygon approximates a circle. The
    output can be input to a matplotlib.line2D object to create a
    circle to draw around a point on a map.

    Return <xvals>, <yvals>.
    """
    delta_theta=2*np.pi/float(npts)
    theta=[xx*delta_theta for xx in range(npts)]
    xvals=[x_centre+radius*np.cos(xx) for xx in theta]
    yvals=[y_centre+radius*np.sin(xx) for xx in theta]
    if verbose:
        print('# create_regular_polygon')
        print('x_centre,y_centre,radius,npts: {0!s}, {1!s}, {2!s}, {3!s}'.format(x_centre,y_centre,radius,npts))
        print('delta_theta: {0!s}'.format(delta_theta))
        print('theta: {0!s}'.format(theta))
        print('xvals: {0!s}'.format(xvals))
        print('yvals: {0!s}'.format(yvals))
    return xvals,yvals

