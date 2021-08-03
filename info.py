"""Information on locations, sections, etc."""

# lon, lat of key locations
locations={}
locations['Jakarta']={'lon':106.7, 'lat':-6.1}
locations['Pameungpeuk']={'lon':107.68, 'lat':-7.65}
locations['Christmas_Island']={'lon':105.6, 'lat':-10.5}
locations['Kota_Kinabalu']={'lon':116.1, 'lat':5.9}
locations['Swallow_Reef']={'lon':113.8, 'lat':7.4}
locations['Darwin']={'lon':130.9, 'lat':-12.4}
locations['Darwin_NW1']={'lon':127.5, 'lat':-11.0}
locations['Darwin_NW2']={'lon':127.0, 'lat':-11.5}
locations['Darwin_W1']={'lon':127.5, 'lat':-12.4}
locations['ELOX1']={'lon':105.9000, 'lat':-10.1667}
locations['ELOX2']={'lon':106.5333, 'lat':-9.2500}

def getlonlat(xx):
    return locations[xx]['lon'],locations[xx]['lat']

# sections (start and end point locations)
sections={}
sections[1]={'location1':'Christmas_Island', 'location2':'Jakarta'}
sections[2]={'location1':'Darwin_W1', 'location2':'Darwin'}
sections[3]={'location1':'Darwin_NW2', 'location2':'Darwin'}

# spatial domains (SDOMAIN; 3 characters)
sdomains={}
sdomains['plp']={'lon1':115, 'lon2':130, 'lat1':5, 'lat2':20, 'desc':'Philippines'}
sdomains['mcw']={'lon1':90, 'lon2':130, 'lat1':-15, 'lat2':20, 'desc':'Western Maritime Continent'}
sdomains['mts']={'lon1':290, 'lon2':305, 'lat1':10, 'lat2':25, 'desc':'Caribbean including Montserrat'}
sdomains['mt2']={'lon1':296, 'lon2':300, 'lat1':15, 'lat2':19, 'desc':'Montserrat'}
sdomains['npl']={'lon1':75, 'lon2':90, 'lat1':20, 'lat2':35, 'desc':'Nepal'}
sdomains['np2']={'lon1':84, 'lon2':88, 'lat1':26, 'lat2':30, 'desc':'Nepal region 2'}
sdomains['bar']={'lon1':302, 'lon2':304, 'lat1':13, 'lat2':15, 'desc':'Eureca4 glider near Barbados'}
sdomains['trm']={'lon1':-180, 'lon2':180, 'lat1':-50, 'lat2':50, 'desc':'TRMM spatial domain'}
sdomains['trp']={'lon1':-180, 'lon2':180, 'lat1':-30, 'lat2':30, 'desc':'Tropics'}

class LevelWidths(object):
    """Calculate level widths from a set of levels. 

    Input is a set of levels, from which half-levels and level widths
    are calculated. Typically the levels correspond to the levels in a
    reanalysis product, or model setup.

    Inputs:

    <levels> a list of levels.

    <level_type> String describing level type. Use same convention as
    in source level_type, e.g., 'plev', 'zlev'.

    <units> e.g. 'hPa'
    """

    def __init__(self,levels,level_type,units):
        """Initialise. """
        self.levels=levels
        self.level_type=level_type
        self.units=units
        self.half_levels_below={}
        self.half_levels_above={}
        self.level_widths={}
        self.nlevels=len(self.levels)
        for index in range(self.nlevels):
            level=self.levels[index]
            if index==0:
                level_below=1e20 # Set to very large value.
                level_above=self.levels[1]
            elif index==self.nlevels-1:
                level_below=self.levels[-2]
                level_above=0
            else:
                level_below=self.levels[index-1]
                level_above=self.levels[index+1]
            half_level_below=0.5*(level_below+level)
            half_level_above=0.5*(level+level_above)
            level_width=half_level_below-half_level_above
            #print(index,level_below,level,level_above,half_level_below,half_level_above,level_width)
            self.level_widths[level]=level_width
            self.half_levels_below[level]=half_level_below
            self.half_levels_above[level]=half_level_above

levels={}
# Level widths for NCEP-DOE reanalysis
levels['ncepdoe']=LevelWidths([1000,925,850,700,600,500,400,300,250,200,150,100,70,50,30,20,10],'plev','hPa')
# Level widths for ERA-Interim reanalysis
levels['erainterim']=LevelWidths([1000,975,950,925,900,875,850,825,800,775,750,700,650,600,550,500,450,400,350,300,250,225,200,175,150,125,100,70,50,30,20,10,7,5,3,2,1],'plev','hPa')

# 'M' type time domains, e.g., 'M0001': parameters
Mtdomain={}
Mtdomain['0001']={'source':'trmm3b42v7p1_sfc_d', 'tseriesfile':'ppt_1_ss_lat_3.125_3.125_lon_101.625_101.625', 'method':'percentile_above', 'percentile':95}
Mtdomain['0002']={'source':'trmm3b42v7p1_sfc_d', 'tseriesfile':'ppt_1_ss_lat_-0.875_-0.875_lon_100.375_100.375', 'method':'percentile_above', 'percentile':95}
Mtdomain['0003']={'source':'trmm3b42v7p1_sfc_d', 'tseriesfile':'ppt_1_ss_lat_-0.875_-0.875_lon_100.375_100.375', 'method':'percentile_above', 'percentile':99}

