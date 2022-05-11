"""Create time domain from index time series."""

import os

import data_analysis as da

#TDOMAINID='rmm006all1'
TDOMAINID='M0003'

#BASEDIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')
BASEDIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')

VERBOSE=True

#==========================================================================

# Create instance of TimeDomain object
#aa=da.TimeDomain(TDOMAINID,verbose=VERBOSE)
aa=da.TimeDomain(TDOMAINID,verbose=VERBOSE,tseriesdir=BASEDIR)

# Create time domain
aa.f_create_time_domain_from_indices()
