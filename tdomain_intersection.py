"""Create new time domain from intersection or other of two other time domains.

Two input time domains
TDOMAINID1: e.g., crossing times of CCKW events at the location
TDOMAINID2: e.g., datetimes of HIW events at a location

Output time domain
TDOMAINID3: common date time from the two input time domains

"""

import copy
import datetime
import os

import cftime

import data_analysis as da

# This is real hack code now. Check output very carefully.

# Choose two time domains. They can be, e.g.,
# time domain for EW crossing times at location
# time domain for dates of extreme precip at location
#TDOMAINID1='CCER103Elat-15-15-ndj98-18-0.05-00UTC' 
TDOMAINID1='CCER103Elat-15-15-ndj98-18-0.05-00UTC-and-vwnd-850-lt-0.5' 
TDOMAINID2='M0004' 

METHOD='and' # 'and', '1not2', '2not1'
LAG=True # Set True to introduce lag offsets
if LAG:
    LAG1=-1
    LAG2=1

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
    # Ad hoc code below. Uncomment the code you need
    #
    if not LAG:
        # exact match of datetimes
        #dt3=[xx for xx in td1.datetimes if xx in td2.datetimes]
        # datetimes just need to be on same day
        dt3=[xx for xx in td1.datetimes if [cftime.DatetimeGregorian(xx[0].year,xx[0].month,xx[0].day)] in td2.datetimes]
        tdomainid3=TDOMAINID1+'-and-'+TDOMAINID2
    else:
        # Include date if a time in TDOMAINID2 occurs within (LAG1,LAG2) of a time in TDOMAINID1
        # Corollary. Lagged composites based on TDOMAINID3 will be centered (i.e., lag 0) on the events
        #   in TDOMAINID2
        dt3=[]
        for xx in td1.datetimes:
            #xx1=xx
            xx1=cftime.DatetimeGregorian(xx[0].year,xx[0].month,xx[0].day)
            for lagc in range(LAG1,LAG2+1):
                xx2=xx1+datetime.timedelta(days=lagc)
                for yy in td2.datetimes:
                    if yy[0]==xx2 and yy not in dt3:
                        dt3.append(yy)
        tdomainid3=TDOMAINID1+'-and-'+TDOMAINID2+'-lag'+str(LAG1)+'-'+str(LAG2)
elif METHOD=='1not2':
    dt3=[xx for xx in td1.datetimes if xx not in td2.datetimes]
    tdomainid3=TDOMAINID2+'-and-not-'+TDOMAINID1
    if LAG:
        raise UserWarning('Not coded this up yet.')
elif METHOD=='2not1':
    if not LAG:
        dt3=[xx for xx in td2.datetimes if xx not in td1.datetimes]
        tdomainid3=TDOMAINID2+'-and-not-'+TDOMAINID1
    else:
        # Set dt3 initially to be a copy of td2.datetimes, then loop over datetimes
        # in dt1 (with lags) and remove datetime from dt3 if there is overlap.
        dt3=copy.copy(td2.datetimes)
        for xx in td1.datetimes:
            #xx1=xx
            xx1=cftime.DatetimeGregorian(xx[0].year,xx[0].month,xx[0].day)
            for lagc in range(LAG1,LAG2+1):
                xx2=xx1+datetime.timedelta(days=lagc)
                for yy in dt3:
                    if yy[0]==xx2:
                        dt3.remove(yy)
        tdomainid3=TDOMAINID2+'-and-not-'+TDOMAINID1+'-lag'+str(LAG1)+'-'+str(LAG2)
else:
    raise UserWarning('Invalid METHOD.')

# Create time domain 3
td3=da.TimeDomain(tdomainid3,verbose=VERBOSE)
td3.datetimes=dt3
td3.datetime2ascii()
header0=td1.header[0][:-1]+' : '+td2.header[0][1:-1]+' : '+METHOD+'\n'
if not LAG:
    header1=td1.header[1][:-1]+' : '+td2.header[1][1:]
else:
    header1=td1.header[1][:-1]+' : '+td2.header[1][1:-1]+' : lag '+str(LAG1)+' to '+str(LAG2)
td3.header=[header0,header1]
td3.write_ascii()

# Check. Count number of events in time domains
td1.f_nevents()
td2.f_nevents()
td3.f_nevents()
print('td1.nevents: {0.nevents!s}'.format(td1))
print('td2.nevents: {0.nevents!s}'.format(td2))
print('td3.nevents: {0.nevents!s}'.format(td3))
