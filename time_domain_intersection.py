"""Create new time domain from intersection or other of two other time domains.

Two input time domains
TDOMAINID1: e.g., crossing times of CCKW events at the location
TDOMAINID2: e.g., datetimes of HIW events at a location

Output time domain
TDOMAINID3: common date time from the two input time domains

"""

import os

import datetime

import data_analysis as da

TDOMAINID1='CCEK102E98-20-3h' # time domain for CCKW crossing times at location
TDOMAINID2='M0002' # time domain for dates of extreme precip at location

METHOD='and' # 'and', '1not2', '2not1'

VERBOSE=True

#--------------------------------------------------------------------------------------

# Create time domain objects
td1=da.TimeDomain(TDOMAINID1,verbose=VERBOSE)
td2=da.TimeDomain(TDOMAINID2,verbose=VERBOSE)

# Read in time domains
td1.read_ascii()
td1.ascii2datetime()
td2.read_ascii()
td2.ascii2datetime()

# Create list of intersecting datetimes
if METHOD=='and':
    # exact match of datetimes
    #dt3=[xx for xx in td1.datetimes if xx in td2.datetimes]
    # datetimes just need to be on same day
    dt3=[xx for xx in td1.datetimes if [datetime.datetime(xx[0].year,xx[0].month,xx[0].day)] in td2.datetimes]
    tdomainid3=TDOMAINID1+'-and-'+TDOMAINID2
elif METHOD=='1not2':
    dt3=[xx for xx in td1.datetimes if xx not in td2.datetimes]
    tdomainid3=TDOMAINID2+'-and-not-'+TDOMAINID1
elif METHOD=='2not1':
    dt3=[xx for xx in td2.datetimes if xx not in td1.datetimes]
    tdomainid3=TDOMAINID2+'-and-not-'+TDOMAINID1
else:
    raise UserWarning('Invalid METHOD.')

# Create time domain 3
td3=da.TimeDomain(tdomainid3,verbose=VERBOSE)
td3.datetimes=dt3
td3.datetime2ascii()
header0=td1.header[0][:-1]+' : '+td2.header[0][1:-1]+' : '+METHOD+'\n'
header1=td1.header[1][:-1]+' : '+td2.header[1][1:]
td3.header=[header0,header1]
td3.write_ascii()

td1.f_nevents()
td2.f_nevents()
td3.f_nevents()
print('td1.nevents: {0.nevents!s}'.format(td1))
print('td2.nevents: {0.nevents!s}'.format(td2))
print('td3.nevents: {0.nevents!s}'.format(td3))
