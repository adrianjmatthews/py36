"""Create time domain from index time series."""

import os

import data_analysis as da

TDOMAINID='rmm008djf1'
#TDOMAINID='M0005'

VERBOSE=True

#==========================================================================

# Create instance of TimeDomain object
aa=da.TimeDomain(TDOMAINID,verbose=VERBOSE)

# Create time domain
aa.f_create_time_domain_from_indices()
