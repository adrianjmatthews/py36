"""Combine latitudes into single cube using data_analysis.CombineLatitudes.

Best to run as a single script through job_ada.sub. Takes approx. 1 day.

Do not be tempted to loop from YEAR_BEG to YEAR_END in run_scripts_sub.py
Each year run would have to read in the same TIME1-TIME2 full length files and it
slows the whole system down.
"""

import os

import cftime
import datetime
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt

import data_analysis as da
import info

BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')

ARCHIVE=True
BASEDIR_ARCHIVE=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

VAR_NAME='uwnd'; 
LEVEL=850; SOURCE1='era5gloerai_plev_3h';  SOURCE2='era5gloeraiER1_plev_3h'

LATITUDES=['-89.4629', '-88.7669', '-88.0669', '-87.3660', '-86.6648', '-85.9633', '-85.2618', '-84.5602', '-83.8586', '-83.1569', '-82.4553', '-81.7536', '-81.0519', '-80.3502', '-79.6485', '-78.9468', '-78.2450', '-77.5433', '-76.8416', '-76.1399', '-75.4381', '-74.7364', '-74.0346', '-73.3329', '-72.6312', '-71.9294', '-71.2277', '-70.5260', '-69.8242', '-69.1225', '-68.4207', '-67.7190', '-67.0172', '-66.3155', '-65.6137', '-64.9120', '-64.2102', '-63.5085', '-62.8067', '-62.1050', '-61.4033', '-60.7015', '-59.9998', '-59.2980', '-58.5963', '-57.8945', '-57.1928', '-56.4910', '-55.7893', '-55.0875', '-54.3858', '-53.6840', '-52.9823', '-52.2805', '-51.5788', '-50.8770', '-50.1753', '-49.4735', '-48.7718', '-48.0700', '-47.3683', '-46.6665', '-45.9647', '-45.2630', '-44.5612', '-43.8595', '-43.1577', '-42.4560', '-41.7542', '-41.0525', '-40.3507', '-39.6490', '-38.9472', '-38.2455', '-37.5437', '-36.8420', '-36.1402', '-35.4385', '-34.7367', '-34.0350', '-33.3332', '-32.6315', '-31.9297', '-31.2280', '-30.5262', '-29.8244', '-29.1227', '-28.4209', '-27.7192', '-27.0174', '-26.3157', '-25.6139', '-24.9122', '-24.2104', '-23.5087', '-22.8069', '-22.1052', '-21.4034', '-20.7017', '-19.9999', '-19.2982', '-18.5964', '-17.8947', '-17.1929', '-16.4911', '-15.7894', '-15.0876', '-14.3859', '-13.6841', '-12.9824', '-12.2806', '-11.5789', '-10.8771', '-10.1754', '-09.4736', '-08.7719', '-08.0701', '-07.3684', '-06.6666', '-05.9649', '-05.2631', '-04.5613', '-03.8596', '-03.1578', '-02.4561', '-01.7543', '-01.0526', '-00.3508', '00.3508', '01.0526', '01.7543', '02.4561', '03.1578', '03.8596', '04.5613', '05.2631', '05.9649', '06.6666', '07.3684', '08.0701', '08.7719', '09.4736', '10.1754', '10.8771', '11.5789', '12.2806', '12.9824', '13.6841', '14.3859', '15.0876', '15.7894', '16.4911', '17.1929', '17.8947', '18.5964', '19.2982', '19.9999', '20.7017', '21.4034', '22.1052', '22.8069', '23.5087', '24.2104', '24.9122', '25.6139', '26.3157', '27.0174', '27.7192', '28.4209', '29.1227', '29.8244', '30.5262', '31.2280', '31.9297', '32.6315', '33.3332', '34.0350', '34.7367', '35.4385', '36.1402', '36.8420', '37.5437', '38.2455', '38.9472', '39.6490', '40.3507', '41.0525', '41.7542', '42.4560', '43.1577', '43.8595', '44.5612', '45.2630', '45.9647', '46.6665', '47.3683', '48.0700', '48.7718', '49.4735', '50.1753', '50.8770', '51.5788', '52.2805', '52.9823', '53.6840', '54.3858', '55.0875', '55.7893', '56.4910', '57.1928', '57.8945', '58.5963', '59.2980', '59.9998', '60.7015', '61.4033', '62.1050', '62.8067', '63.5085', '64.2102', '64.9120', '65.6137', '66.3155', '67.0172', '67.7190', '68.4207', '69.1225', '69.8242', '70.5260', '71.2277', '71.9294', '72.6312', '73.3329', '74.0346', '74.7364', '75.4381', '76.1399', '76.8416', '77.5433', '78.2450', '78.9468', '79.6485', '80.3502', '81.0519', '81.7536', '82.4553', '83.1569', '83.8586', '84.5602', '85.2618', '85.9633', '86.6648', '87.3660', '88.0669', '88.7669', '89.4629']

# Parameters from wheelerkiladis.py that determine input file names
WAVE_TYPE='ER'
TIME1=cftime.DatetimeGregorian(1998,1,1)
TIME2=TIME1+datetime.timedelta(25*365+6-1)-datetime.timedelta(seconds=1)
#TIME1=cftime.DatetimeGregorian(2009,1,1)
#TIME2=TIME1+datetime.timedelta(3*365-1)-datetime.timedelta(seconds=1)

FILEPRE='' # e.g., '', '_rac', '_rac_b20_200_n241', '_rac_rm5_n5'

# Do not be tempted to loop from YEAR_BEG to YEAR_END in run_scripts_sub.py
# Each year run would have to read in the same TIME1-TIME2 full length files and it
# slows the whole system down.
YEAR_BEG=2021; YEAR_END=2022

PLOT=False

VERBOSE=2

#------------------------------------------------------------------

descriptor={}
descriptor['verbose']=VERBOSE
descriptor['basedir']=BASEDIR
descriptor['archive']=ARCHIVE
descriptor['basedir_archive']=BASEDIR_ARCHIVE
descriptor['source1']=SOURCE1
descriptor['source2']=SOURCE2
descriptor['var_name']=VAR_NAME
descriptor['level']=LEVEL
descriptor['latitudes']=LATITUDES
descriptor['filepre']=FILEPRE
descriptor['wave_type']=WAVE_TYPE
descriptor['time1']=TIME1
descriptor['time2']=TIME2

# Create instance of CombineLatitudes object
aa=da.CombineLatitudes(**descriptor)

for year in range(YEAR_BEG,YEAR_END+1):
    print('### year={0!s}'.format(year))
    aa.year=year
    aa.f_combine_latitudes()

if PLOT:
    fig=plt.figure()
    x1=aa.data_all
    x2=x1.extract(iris.Constraint(latitude=0.0))
    qplt.contourf(x2,yrev=1)
    
    plt.show()
    fig.savefig('/gpfs/home/e058/tmp/fig1.png')
