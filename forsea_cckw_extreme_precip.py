"""Ad hoc script to calculate dependence of extreme precip at location on CCKW existence."""

import datetime

import data_analysis as da

TDOMAINID1='CCEK102E98-20-00UTC' # time domain for CCKW crossing times at location
TDOMAINID2='M0003' # time domain for dates of extreme precip at location

# (Start,end) valid time range for subsequent analysis
TIME1=datetime.datetime(1998,1,1)
TIME2=datetime.datetime(2020,12,31)

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

# Extract datetimes from time domain in valid time range
def extract_datetimes(td,time1,time2):
    return [xx[0] for xx in td.datetimes if time1<=xx[0]<=time2]
dtv1=extract_datetimes(td1,TIME1,TIME2)
dtv2=extract_datetimes(td2,TIME1,TIME2)

# Calculate total number of days in valid time range
x1=TIME2-TIME1
ndays=x1.days
print('ndays: {0!s}'.format(ndays))

# Calculate statistics
ndays1=len(dtv1)
ndays2=len(dtv2)
print('ndays1: {0!s}'.format(ndays1))
print('ndays2: {0!s}'.format(ndays2))
#
def overlapping_days(dt1,dt2,nday):
    days_overlap=[]
    timedelta1=datetime.timedelta(days=nday)
    for dtc in dt1:
        overlap=False
        for xx in dt2:
            if xx-timedelta1 <= dtc <= xx+timedelta1:
                overlap=True
        if overlap:
            days_overlap.append(dtc)
    return days_overlap
overlap0=overlapping_days(dtv1,dtv2,0)
overlap1=overlapping_days(dtv1,dtv2,1)
overlap2=overlapping_days(dtv1,dtv2,2)
overlap3=overlapping_days(dtv1,dtv2,3)
noverlap0=len(overlap0)
noverlap1=len(overlap1)
noverlap2=len(overlap2)
noverlap3=len(overlap3)
print('noverlap0: {0!s}'.format(noverlap0))
print('noverlap1: {0!s}'.format(noverlap1))
print('noverlap2: {0!s}'.format(noverlap2))
print('noverlap3: {0!s}'.format(noverlap3))
#
print('Define event A as observing extreme precipitation on a given day.')
print('Define event B as observing a CCKW  on a given day.')
pA=ndays2/ndays
print('Unconditional probability of observing extreme rainfall, p(A): {0!s}'.format(pA))
pB=ndays1/ndays
print('Unconditional probability of observing CCKW on a given day, p(B): {0!s}'.format(pB))
pAandB=noverlap0/ndays
print('Probability of observing extreme rainfall and CCKW on same day, p(A^B): {0!s}'.format(pAandB))
pAgivenB=pAandB/pB
print('Conditional probability of observing extreme rainfall given there is a CCKW, p(A/B): {0!s}'.format(pAgivenB))
#
days_CCKWplusminus1=[]
dtc=TIME1
timedelta1=datetime.timedelta(days=1)
while dtc<=TIME2:
    include=False
    for xx in dtv1:
        if xx-timedelta1 <= dtc <= xx+timedelta1:
            include=True
    if include:
        days_CCKWplusminus1.append(dtc)
    dtc=dtc+timedelta1
ndays_CCKWplusminus1=len(days_CCKWplusminus1)
print('ndays_CCKWplusminus1: {0!s}'.format(ndays_CCKWplusminus1))
print('Define event C as observing a CCKW within plus or minus 1 day, on a given day.')
pC=ndays_CCKWplusminus1/ndays
print('Unconditional probability of observing CCKW within plus or minus 1 day on a given day, p(C): {0!s}'.format(pC))
pAandC=noverlap1/ndays
print('Probability of observing extreme rainfall and being within plus or minus 1 day of CCKW, p(A^C): {0!s}'.format(pAandC))
pAgivenC=pAandC/pC
print('Conditional probability of observing extreme rainfall given there is a CCKW within plus or minus 1 day, p(A/C): {0!s}'.format(pAgivenC))
inflationfactor=pAgivenC/pA
print('Inflationary factor: p(A/C) / p(A): {0!s}'.format(inflationfactor))
