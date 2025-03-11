"""Submits python scripts in parallel on ada using SLURM.

Takes the python script in FILENAME and runs it many times, in
parallel, by nested loops over LOOPVAR1, LOOPVAR2, LOOPVAR3, LOOPVAR4,
LOOPVAR5, which represent e.g., YEAR, MONTH, TDOMAINID, VAR_NAME,
LEVEL

If one of the LOOPVAR variables is not needed it is set to a dummy
list of [-999]

Make sure that any 'YEAR=...', 'LEVEL=...' assignments in FILENAME are
commented out.

Also, make sure that PLOT is set to False.

For each instance, i.e., each value of YEAR, LEVEL, etc, a temporary
directory is created.  A copy of FILENAME with the relevant 'YEAR=...'
and 'LEVEL=...' statements at the beginning is made in the temporary
directory.  Then a 'sub' file is made in the temporary directory, and
the sub file is submitted as a batch job.

'squeue -u <userid>' to look at jobs, once submitted.

'scancel 0' to kill all jobs (CHECK THIS WORKS).

sacct -j <jobid> --format='JobID,MaxRSS,reqmem,CPUTime,Timelimit,elapsed' # to check resources used by job
sacct -e # get list of all fields to get values for

If note sure how much memory a job will take, set memory to a large
value, run one job, check memory usage using sacct, then set memory
variable below for future jobs. Running a job with not enough memory
requested will cause swap space to be used, which has a highly
detrimental effect on performance and time taken.

"""

import os
import shutil
import time

#FILENAME='test.py'
#FILENAME='preprocess.py'
#FILENAME='zero_missing_values.py'
FILENAME='anncycle.py'
#FILENAME='regrid.py'
#FILENAME='time_average.py'
#FILENAME='filter.py'
#FILENAME='spatial_subset.py'
#FILENAME='mean.py'
#FILENAME='lagged_mean.py'
#FILENAME='vwndptap.py'
#FILENAME='vrtbudget.py'
#FILENAME='vrtbudget_combine.py'
#FILENAME='truncate.py'
#FILENAME='tdomain_create.py'
#FILENAME='diurnal_cycle.py'
#FILENAME='imerg_get.py'
#FILENAME='imerg_hdf5tonetcdf.py'
#FILENAME='imerg2trmm.py'
#FILENAME='nc2txt_time_series.py'
#FILENAME='wind.py'
#FILENAME='overturning_potential.py'
#FILENAME='overturning_vector_streamfunction.py'
#FILENAME='omega_decomposition.py'
#FILENAME='mass_flux_decomposition.py'
#FILENAME='wheelerkiladis.py'
#FILENAME='combine_latitudes.py'
#FILENAME='tmp_reorder.py'
#FILENAME='ostia_get.py'
#FILENAME='wndspd.py'
#FILENAME='tdomain_from_lagrangian.py'
#FILENAME='copernicus_get.py'
#FILENAME='subtract.py'
#FILENAME='remove_coord.py'
#FILENAME='plot_topography.py'

# Experiment name for temporary sub directories
TMPEXP='tmp0_'+FILENAME.split('.')[0] 
DIR1=os.path.join(os.path.sep,'gpfs','home','e058','tmp','sub_'+TMPEXP)

NICE=0 # Set to 10 to lower priority so new interactive sessions can be
        # started. If ada becomes full and cannot get batch jobs to run, set to 0.

# Set whether or not to run each script with exclusive flag.
# Only use if having issues running out of memory
# If true, will use all memory on node, but there is a max of ~8(?) exclusive 
#   jobs allowed per user, so do not set for large loops unless absolutely necessary
# Think this is equivalent to requesting 64G memory (on compute-16-64 queue)
# NB with exclusive flag, also use --mem=0 to allocate all memory
EXCLUSIVE=False
if FILENAME in ['combine_latitudes.py']:
    #EXCLUSIVE=True # Use all memory on node
    pass

# Set number of maximum jobs that can run simultaneously (set to 0 to disable).
# Can be useful if system gets overwhelmed with i/o with many jobs running.
# This is used with the --dependency=singleton option for sbatch later.
# NB If this still causes problems consider using slurm job array submission.
NMAXJOBS=0

# Set queue, memory, runtime requirements
# Default: RAM per core on compute-16-64 queue is 4 Gb. The 16 in compute-16-64
# is number of nodes available. Hence 16x4=64 GB max memory. 
# In reality, memory request should be set to slightly below this, e.g., 60G.
# queues are compute-16-64, compute-24-96, compute-24-128
# Set default values then override below if required
#  (dedicated high throughput queue)
queue='compute-64-512,compute-24-96' # new, use these
memory='12G'
runtime='6:00:00'
if FILENAME in ['mean.py','wndspd.py','nc2txt_time_series.py']:
    # mean.py, lagged_mean.py (with lazy_load True) ~1G
    #memory='4G'
    # spatial_subset.py, wheelerkiladis.py ~0.7G but timed out at 4 hr with only 4G requested
    #   yet completed in 5 min with 12G requested.
    # Recommend using minimum of 12G; doesn't seem to delay job running.
    pass
elif FILENAME in ['filter.py']:
    # filter.py ~9G erainterim_plev_d with splitblock true
    memory='70G' 
elif FILENAME in ['lagged_mean.py','tdomain_create.py']:
    runtime='24:00:00'
elif FILENAME in ['wheelerkiladis.py']:
    # if it doesn't run in 30 minutes it won't run. i/o issue.
    # used 10G of MaxRSS, but failed with memory of 12G!
    memory='24G' 
elif FILENAME in ['anncycle.py']:
    memory='144G' 
elif FILENAME in ['vrtbudget.py']:
    # vrtbudget.py seems to accumulate memory when looped over years within
    # its own script and fails with out of memory.
    # Run it here by looping over years
    memory='42G'
    runtime='24:00:00'
elif FILENAME in ['combine_latitudes.py']:
    memory='80G'
    runtime='36:00:00'
elif FILENAME in ['omega_decomposition.py']:
    memory='52G'
elif FILENAME in ['preprocess.py','time_average.py']:
    # preprocess.py 100G for imergtrm_sfc_30m.
    # time_average.py 27G for imergtrm_sfc_30m.
    queue='compute-64-512'
    memory='100G' # imergtrm
    #memory='24G'
    runtime='18:00:00'
elif FILENAME in ['regrid.py']:
    # regrid.py 48G for era5gloerai
    #queue='compute-64-512' # Use this queue for imerg preprocess
    #memory='100G' # imergtrm
    memory='100G'
    runtime='18:00:00'
elif FILENAME in ['copernicus_get.py','imerg_get.py']:
    memory='48G'
    runtime='72:00:00'
elif FILENAME in ['plot_topography.py']:
    memory='48G'

# 'var1': 'YEAR' or 'YEAR_BEG'
# 'var2': 'MONTH' or 'SOURCE'
# 'var3': 'TDOMAINID' or 'BAND1_VAL1' or 'LAT1'
# 'var4': 'VAR_NAME' or 'BAND2_VAL1'
# 'var5': 'LEVEL' or 'LONC'
VARDICT={'var1':'XXX',
         'var2':'MONTH',
         'var3':'XXX',
         'var4':'VAR_NAME',
         'var5':'LEVEL'}
if FILENAME in ['copernicus_get.py','combine_latitudes.py']:
    VARDICT['var1']='YEAR_BEG'
else:
    VARDICT['var1']='YEAR'
if FILENAME=='spatial_subset.py':
    VARDICT['var3']='BAND1_VAL1'
elif FILENAME=='wheelerkiladis.py':
    VARDICT['var3']='LAT1'
else:
    VARDICT['var3']='TDOMAINID'

# Initial dummy values, to be overwritten if needed
LOOPVAR1=LOOPVAR2=LOOPVAR3=LOOPVAR4=LOOPVAR5=['X']

######################################################
#LOOPVAR1=range(2003,2020+1)
#LOOPVAR1=range(1998,2022+1,9) # copernicus_get.py
#LOOPVAR1=range(2001,2023+1,8) # copernicus_get.py for Natasha YEAREND=YEARBEG+7
#LOOPVAR1=range(1998,2022+1,5) # combine_latitudes.py

######################################################
#LOOPVAR2=range(1,12+1) # Set MONTH ranges if outfile_frequency is less than 'year'
#LOOPVAR2=[3,]
#LOOPVAR2=['sg579m031oi01_zlev_h','sg534m031oi01_zlev_h','sg532m031oi01_zlev_h','sg620m031oi01_zlev_h','sg613m031oi01_zlev_h']
#LOOPVAR2=range(100) # NODE_NUMBER: Number of processors to farm out mean.py for Monte Carlo simulation of null distribution

######################################################
#LOOPVAR3=['jan01-22','feb01-22','mar01-22','apr01-22','may01-22','jun01-22','jul01-22','aug01-22','sep01-22','oct01-22','nov01-22','dec01-22']
#LOOPVAR3=['ann'+str(xx).zfill(4) for xx in range(1979,2017+1)]
#LOOPVAR3=['rmm007n2a'+str(xx).zfill(1) for xx in range(1,8+1) ] + ['rmm007m2o'+str(xx).zfill(1) for xx in range(1,8+1) ]
#LOOPVAR3=['djf0001-2223','mam01-23','jja01-23','son00-22']
#LOOPVAR3=['CCER103Elat-15-15-ndj98-18-0.05-00UTC-and-vwnd-850-lt-0.5','CCER103Elat-15-15-ndj98-18-0.05-00UTC-and-vwnd-850-lt-0.5-and-M0004-lag-1-1']
#LOOPVAR3=['CCEK-30E98-20-0.3-00UTC','CCEK-15E98-20-0.3-00UTC','CCEK0E98-20-0.3-00UTC','CCEK30E98-20-0.3-00UTC']
#LOOPVAR3=['CCEK102E98-18-00UTC-and-M0002b','M0002b-and-not-CCEK102E98-18-00UTC'] # Natasha
#LOOPVAR3=['CCEK'+str(xx)+'E98-20-0.3-00UTC' for xx in range(-170,180,10)]
#LOOPVAR3=['01.0526']
#LOOPVAR3=['-89.4629', '-88.7669', '-88.0669', '-87.3660', '-86.6648', '-85.9633', '-85.2618', '-84.5602', '-83.8586', '-83.1569', '-82.4553', '-81.7536', '-81.0519', '-80.3502', '-79.6485', '-78.9468', '-78.2450', '-77.5433', '-76.8416', '-76.1399', '-75.4381', '-74.7364', '-74.0346', '-73.3329', '-72.6312', '-71.9294', '-71.2277', '-70.5260', '-69.8242', '-69.1225', '-68.4207', '-67.7190', '-67.0172', '-66.3155', '-65.6137', '-64.9120', '-64.2102', '-63.5085', '-62.8067', '-62.1050', '-61.4033', '-60.7015', '-59.9998', '-59.2980', '-58.5963', '-57.8945', '-57.1928', '-56.4910', '-55.7893', '-55.0875', '-54.3858', '-53.6840', '-52.9823', '-52.2805', '-51.5788', '-50.8770', '-50.1753', '-49.4735', '-48.7718', '-48.0700', '-47.3683', '-46.6665', '-45.9647', '-45.2630', '-44.5612', '-43.8595', '-43.1577', '-42.4560', '-41.7542', '-41.0525', '-40.3507', '-39.6490', '-38.9472', '-38.2455', '-37.5437', '-36.8420', '-36.1402', '-35.4385', '-34.7367', '-34.0350', '-33.3332', '-32.6315', '-31.9297', '-31.2280', '-30.5262', '-29.8244', '-29.1227', '-28.4209', '-27.7192', '-27.0174', '-26.3157', '-25.6139', '-24.9122', '-24.2104', '-23.5087', '-22.8069', '-22.1052', '-21.4034', '-20.7017', '-19.9999', '-19.2982', '-18.5964', '-17.8947', '-17.1929', '-16.4911', '-15.7894', '-15.0876', '-14.3859', '-13.6841', '-12.9824', '-12.2806', '-11.5789', '-10.8771', '-10.1754', '-09.4736', '-08.7719', '-08.0701', '-07.3684', '-06.6666', '-05.9649', '-05.2631', '-04.5613', '-03.8596', '-03.1578', '-02.4561', '-01.7543', '-01.0526', '-00.3508', '00.3508', '01.0526', '01.7543', '02.4561', '03.1578', '03.8596', '04.5613', '05.2631', '05.9649', '06.6666', '07.3684', '08.0701', '08.7719', '09.4736', '10.1754', '10.8771', '11.5789', '12.2806', '12.9824', '13.6841', '14.3859', '15.0876', '15.7894', '16.4911', '17.1929', '17.8947', '18.5964', '19.2982', '19.9999', '20.7017', '21.4034', '22.1052', '22.8069', '23.5087', '24.2104', '24.9122', '25.6139', '26.3157', '27.0174', '27.7192', '28.4209', '29.1227', '29.8244', '30.5262', '31.2280', '31.9297', '32.6315', '33.3332', '34.0350', '34.7367', '35.4385', '36.1402', '36.8420', '37.5437', '38.2455', '38.9472', '39.6490', '40.3507', '41.0525', '41.7542', '42.4560', '43.1577', '43.8595', '44.5612', '45.2630', '45.9647', '46.6665', '47.3683', '48.0700', '48.7718', '49.4735', '50.1753', '50.8770', '51.5788', '52.2805', '52.9823', '53.6840', '54.3858', '55.0875', '55.7893', '56.4910', '57.1928', '57.8945', '58.5963', '59.2980', '59.9998', '60.7015', '61.4033', '62.1050', '62.8067', '63.5085', '64.2102', '64.9120', '65.6137', '66.3155', '67.0172', '67.7190', '68.4207', '69.1225', '69.8242', '70.5260', '71.2277', '71.9294', '72.6312', '73.3329', '74.0346', '74.7364', '75.4381', '76.1399', '76.8416', '77.5433', '78.2450', '78.9468', '79.6485', '80.3502', '81.0519', '81.7536', '82.4553', '83.1569', '83.8586', '84.5602', '85.2618', '85.9633', '86.6648', '87.3660', '88.0669', '88.7669', '89.4629']
LOOPVAR3=['87.8638', '85.0965', '82.3129', '79.5256', '76.7369', '73.9475', '71.1578', '68.3678', '65.5776', '62.7874', '59.9970', '57.2066', '54.4162', '51.6257', '48.8352', '46.0447', '43.2542', '40.4637', '37.6731', '34.8825', '32.0920', '29.3014', '26.5108', '23.7202', '20.9296', '18.1390', '15.3484', '12.5578', '9.7671', '6.9765', '4.1859', '1.3953', '-1.3953', '-4.1859', '-6.9765', '-9.7671', '-12.5578', '-15.3484', '-18.1390', '-20.9296', '-23.7202', '-26.5108', '-29.3014', '-32.0920', '-34.8825', '-37.6731', '-40.4637', '-43.2542', '-46.0447', '-48.8352', '-51.6257', '-54.4162', '-57.2066', '-59.9970', '-62.7874', '-65.5776', '-68.3678', '-71.1578', '-73.9475', '-76.7369', '-79.5256', '-82.3129', '-85.0965', '-87.8638']
#LOOPVAR3=['16.95','17.05','17.15','17.25','17.35','17.45','17.55']

######################################################
#LOOPVAR4=['uwnd','vwnd','vrt','div']
#LOOPVAR4=['uwnd','vwnd','psfc']
#LOOPVAR4=['dvrtdt','m_uwnd_dvrtdx','m_vwnd_dvrtdy','m_omega_dvrtdp','m_vrt_div','m_ff_div','m_beta_vwnd','m_domegadx_dvwnddp','domegady_duwnddp','source_dvrtdt','res_dvrtdt']
#LOOPVAR4=['vrt_horiz_adv','vrt_stretch','vrt_tilt']
#LOOPVAR4=['m_vwndbar_dvrtdyprime','m_vwndprime_dvrtdybar','m_vwndprime_dvrtdyprime']
#LOOPVAR4=['m_uwndbar_dvrtdxbar','m_uwndbar_dvrtdxprime','m_uwndprime_dvrtdxbar','m_uwndprime_dvrtdxprime','m_vwndbar_dvrtdybar','m_vwndbar_dvrtdyprime','m_vwndprime_dvrtdybar','m_vwndprime_dvrtdyprime']
#LOOPVAR4=['lat','lon','tsc','sa']
#LOOPVAR4=['121.75','121.85','121.95','122.05','122.15','122.25','122.35']

######################################################
#LOOPVAR5=[1000,975,950,925,900,875,850,825,800,775,750,700,650,600,550,500,450,400,350,300,250,225,200,175,150,125,100] # era5
#LOOPVAR5=[1000,925,850,700,600,500,400,300,250,200,150,100,70,50,30,20,10]
#LOOPVAR5=[200,500,850,975]
#LOOPVAR5=[-80]
#LOOPVAR5=[xx for xx in range(60,180,5)]; LOOPVAR5.remove(75)
#LOOPVAR5=[xx for xx in range(-100,180,10)] # CCKW longitudes

#--------------------------------------------------------------------

print('# Check all assignments of LOOPVAR variables are commented out in script')
def check_commented_out(loopvar,var,lines):
    print(loopvar,var)
    if loopvar!=['X']:
        for linec in lines:
            if linec[0]!='#' and var+'=' in linec and var+'==' not in linec and '#'+var+'=' not in linec:
                raise ValueError(var+' is not commented out')
            if 'PLOT=True' in linec:
                raise ValueError('PLOT is set to True')

f1=open(os.path.join(os.path.sep,'gpfs','home','e058','home','PythonScripts','py36',FILENAME))
lines=f1.readlines()
check_commented_out(LOOPVAR1,VARDICT['var1'],lines)
check_commented_out(LOOPVAR2,VARDICT['var2'],lines)
check_commented_out(LOOPVAR3,VARDICT['var3'],lines)
check_commented_out(LOOPVAR4,VARDICT['var4'],lines)
check_commented_out(LOOPVAR5,VARDICT['var5'],lines)

print('# Get initial working directory')
cwd=os.getcwd()

print('# Remove old version of DIR1: {0!s}'.format(DIR1))
try:
    shutil.rmtree(DIR1)
    print('Removed DIR1')
except:
    print('DIR1 does not exist or cannot be removed')

print("# Create DIR1")
os.mkdir(DIR1)

print("# (Nested) Loop over YEAR, LEVEL, VAR etc.")
index=0
scriptfile_original=os.path.join(os.path.sep,'gpfs','home','e058','home','PythonScripts','py36',FILENAME)
njobs=len(LOOPVAR1)*len(LOOPVAR2)*len(LOOPVAR3)*len(LOOPVAR4)*len(LOOPVAR5)
print('njobs: {0!s}'.format(njobs))
for var1 in LOOPVAR1: 
    for var2 in LOOPVAR2: 
        for var3 in LOOPVAR3: 
            for var4 in LOOPVAR4:
                for var5 in LOOPVAR5:
                    print("### index:",index)
                    print('var1: {0!s}'.format(var1))
                    print('var2: {0!s}'.format(var2))
                    print('var3: {0!s}'.format(var3))
                    print('var4: {0!s}'.format(var4))
                    print('var5: {0!s}'.format(var5))
                    dircname=str(var1)+'_'+str(var2)+'_'+str(var3)+'_'+str(var4)+'_'+str(var5)

                    # Create temporary local subdirectory
                    dir2=os.path.join(DIR1,'tmp_'+dircname)
                    os.mkdir(dir2)
                    print("dir2:",dir2)
                    
                    # Copy data_analysis.py etc. to local subdirectory
                    dir99=os.path.join(os.path.sep,'gpfs','home','e058','home','PythonScripts','py36')
                    shutil.copy(os.path.join(dir99,'data_analysis.py'),dir2)
                    shutil.copy(os.path.join(dir99,'mypaths.py'),dir2)
                    shutil.copy(os.path.join(dir99,'plotter.py'),dir2)
                    shutil.copy(os.path.join(dir99,'info.py'),dir2)
            
                    # Create a composite python script that will be run under sbatch
                    # First line(s) sets YEAR, LEVEL, VAR etc.
                    scriptfile_local=os.path.join(dir2,'local_script.py')
                    fout=open(scriptfile_local,'w')
                    # var1
                    key=VARDICT['var1']
                    if type(var1)==type('hello'):
                        fout.write(key+'="'+str(var1)+'"\n')
                    else:
                        fout.write(key+'='+str(var1)+'\n')
                    # var2
                    key=VARDICT['var2']
                    if type(var2)==type('hello'):
                        fout.write(key+'="'+str(var2)+'"\n')
                    else:
                        fout.write(key+'='+str(var2)+'\n')
                    # var3
                    key=VARDICT['var3']
                    if type(var3)==type('hello'):
                        fout.write(key+'="'+str(var3)+'"\n')
                    else:
                        fout.write(key+'='+str(var3)+'\n')
                    # var4
                    key=VARDICT['var4']
                    if type(var4)==type('hello'):
                        fout.write(key+'="'+str(var4)+'"\n')
                    else:
                        fout.write(key+'='+str(var4)+'\n')
                    # var5
                    key=VARDICT['var5']
                    if type(var5)==type('hello'):
                        fout.write(key+'="'+str(var5)+'"\n')
                    else:
                        fout.write(key+'='+str(var5)+'\n')
                    fout.close()
                    # Remainder of script is copy of FILENAME
                    # Make sure that all assignments of YEAR, LEVEL, VAR etc
                    # are commented out in FILENAME
                    os.system('cat '+scriptfile_original+' >> '+scriptfile_local)

                    # Write a local sub file for hpc
                    fout=open(os.path.join(dir2,'subfile'),'w')
                    fout.write('#!/bin/bash\n')
                    fout.write('#SBATCH -t '+runtime+'\n')
                    fout.write('#SBATCH -p '+queue+'\n')
                    fout.write('#SBATCH --mem '+memory+'\n')
                    if NMAXJOBS==0:
                        fout.write('#SBATCH --job-name='+FILENAME[:3]+'_'+str(index)+'\n')
                    else:
                        fout.write('#SBATCH --job-name='+FILENAME[:3]+'_'+str( divmod(index,NMAXJOBS)[1] )+'\n')
                    fout.write('#SBATCH --nice='+str(NICE)+'\n')
                    if EXCLUSIVE:
                        fout.write('#SBATCH --exclusive\n')
                    fout.write('#SBATCH -o out_'+str(index)+'.txt\n')
                    fout.write('#SBATCH -e err_'+str(index)+'.txt\n')
                    if index==njobs-1:
                        # Send email when last job is completed or fails
                        fout.write('#SBATCH --mail-type=END,FAIL\n')
                    else:
                        fout.write('#SBATCH --mail-type=NONE\n')
                    #fout.write('#SBATCH --mail-user=e058@uea.ac.uk\n')
                    fout.write('#SBATCH --mail-user=adrian.j.matthews@gmail.com\n')
                    #
                    #fout.write('#SBATCH --ntasks=1\n')
                    #fout.write('#SBATCH --nodes=1\n')
                    #fout.write('#SBATCH --cpus-per-task=8\n')
                    #
                    # ib nodes
                    #fout.write('#SBATCH --qos=ib\n') # On 2023-05-17, got 'Invalid qos specification'  when trying to submit to ib, so removed this from submission job
                    #
                    #fout.write('#SBATCH --qos=adrian\n')
                    #
                    fout.write('#SBATCH --nodes=1\n')
                    fout.write('#SBATCH --ntasks=1\n')
                    fout.write('#SBATCH --cpus-per-task=1\n')
                    #
                    fout.write('\n')
                    #
                    fout.write('export OMP_NUM_THREADS=1\n')
                    #
                    # If the module load and source activate lines are run, fails with cannot find iris
                    # It seems to be accessing the right set up without these lines. I don't understand how.
                    # Leave them commented out.
                    #fout.write('module load python/anaconda/2019.10/3.7\n')
                    #fout.write('source activate py36\n')
                    #
                    fout.write('python ./local_script.py\n')
                    fout.close()
            
                    # Change working directory to local directory
                    os.chdir(dir2)
                    
                    # Submit or run script
                    if NMAXJOBS!=0:
                        os.system('sbatch --dependency=singleton ./subfile')
                    else:
                        os.system('sbatch ./subfile')
                    
                    # Increment counter
                    index+=1

                    # Pause for a short time
                    #time.sleep(3)

print('# Change back to initial working directory')
os.chdir(cwd)
