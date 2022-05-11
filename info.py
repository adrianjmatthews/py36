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
locations['Padang']={'lon':100.4172, 'lat':-0.9471}

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
sdomains['ewa']={'lon1':-20, 'lon2':30, 'lat1':0, 'lat2':60, 'desc':'Europe and West Africa'}

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

# ERA-Interim hybrid model to pressure level conversion
# Values of A and B parameters from Berrisford et al. (2011)
# Dictionary keys are model levels k (1 to 61)
# Values are a tuple of (A_{k - 1/2}, B_{k - 1/2})
# A values are in hPa, B values are non-dimensional
# Pressure of half levels then calculated as
# p_{k - 1/2} = A_{k - 1/2} + B_{k - 1/2} * p_s
#   where p_s is surface pressure in hPa
# Pressure at full levels is then
# p_k = 0.5* ( p_{k - 1/2} + p_{k + 1/2} )
nlevel_erainterim_model=61
erainterim_abkmalf={ 1: (0.00, 0.00000),
                    2: (0.20, 0.00000),
                    3: (0.38, 0.00000),
                    4: (0.64, 0.00000),
                    5: (0.96, 0.00000),
                    6: (1.34, 0.00000),
                    7: (1.81, 0.00000),
                    8: (2.35, 0.00000),
                    9: (2.98, 0.00000),
                    10: (3.74, 0.00000),
                    11: (4.65, 0.00000),
                    12: (5.76, 0.00000),
                    13: (7.13, 0.00000),
                    14: (8.84, 0.00000),
                    15: (10.95, 0.00000),
                    16: (13.56, 0.00000),
                    17: (16.81, 0.00000),
                    18: (20.82, 0.00000),
                    19: (25.80, 0.00000),
                    20: (31.96, 0.00000),
                    21: (39.60, 0.00000),
                    22: (49.07, 0.00000),
                    23: (60.18, 0.00000),
                    24: (73.07, 0.00000),
                    25: (87.65, 0.00008),
                    26: (103.76, 0.00046),
                    27: (120.77, 0.00182),
                    28: (137.75, 0.00508),
                    29: (153.80, 0.01114),
                    30: (168.19, 0.02068),
                    31: (180.45, 0.03412),
                    32: (190.28, 0.05169),
                    33: (197.55, 0.07353),
                    34: (202.22, 0.09967),
                    35: (204.30, 0.13002),
                    36: (203.84, 0.16438),
                    37: (200.97, 0.20248),
                    38: (195.84, 0.24393),
                    39: (188.65, 0.28832),
                    40: (179.61, 0.33515),
                    41: (168.99, 0.38389),
                    42: (157.06, 0.43396),
                    43: (144.11, 0.48477),
                    44: (130.43, 0.53571),
                    45: (116.33, 0.58617),
                    46: (102.10, 0.63555),
                    47: (88.02, 0.68327),
                    48: (74.38, 0.72879),
                    49: (61.44, 0.77160),
                    50: (49.42, 0.81125),
                    51: (38.51, 0.84737),
                    52: (28.88, 0.87966),
                    53: (20.64, 0.90788),
                    54: (13.86, 0.93194),
                    55: (8.55, 0.95182),
                    56: (4.67, 0.96765),
                    57: (2.10, 0.97966),
                    58: (0.66, 0.98827),
                    59: (0.07, 0.99402),
                    60: (0.00, 0.99763),
                    61: (0.00, 1.00000)}
