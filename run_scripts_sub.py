"""Submits python scripts in parallel on ada using SLURM.

Takes the python script in FILENAME and runs it many times, in
parallel, by nested loops over LOOPVAR1, LOOPVAR2, LOOPVAR3, LOOPVAR4, LOOPVAR5, which represent e.g., YEAR, MONTH, TDOMAINID, VAR_NAME, LEVEL

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

"""

import os
import shutil

#FILENAME='test.py'
#FILENAME='preprocess.py'
#FILENAME='zero_missing_values.py'
#FILENAME='anncycle.py'
#FILENAME='time_average.py'
#FILENAME='filter.py'
#FILENAME='spatial_subset.py'
FILENAME='mean.py'
#FILENAME='time_average.py'
#FILENAME='lagged_mean.py'
#FILENAME='vwndptap.py'
#FILENAME='vrtbudget.py'
#FILENAME='vrtbudget_combine.py'
#FILENAME='truncate.py'
#FILENAME='time_domain_create.py'
#FILENAME='diurnal_cycle.py'
#FILENAME='netcdf2ascii_james.py'
#FILENAME='imerg_wget.py'
#FILENAME='imerg_hdf5tonetcdf.py'
#FILENAME='nc2txt_time_series.py'
#FILENAME='wind.py'
#FILENAME='overturning_potential.py'
#FILENAME='overturning_vector_streamfunction.py'
#FILENAME='omega_decomposition.py'
#FILENAME='mass_flux_decomposition.py'
#FILENAME='wheelerkiladis.py'
#FILENAME='combine_latitudes.py'
#FILENAME='tmp_reorder.py'

# Experiment name for temporary sub directories
TMPEXP='tmp0_'+FILENAME.split('.')[0] 
DIR1=os.path.join(os.path.sep,'gpfs','home','e058','tmp','sub_'+TMPEXP)

NICE=10 # Set to default 10 to lower priority so new interactive sessions can be
        # started. If ada becomes full and cannot get batch jobs to run, set to 0.

# Set whether or not to run each script with exclusive flag.
# Only use if having issues running out of memory
# If true, will use all memory on node, but there is a max of only 5(?) exclusive 
#   jobs allowed per user, so do not set for large loops unless absolutely necessary
if FILENAME in ['combine_latitudes.py']:
    EXCLUSIVE=True # Use all memory on node
else:
    EXCLUSIVE=False

# Set memory requirements
if FILENAME in ['vrtbudget.py']:
    memory='23G'
elif FILENAME in ['omega_decomposition.py']:
    memory='50G'
else:
    memory='5G'

# 'var1': 'YEAR' or 'YEAR_BEG'
# 'var2': 'MONTH' or 'SOURCE'
# 'var3': 'TDOMAINID' or 'BAND1_VAL1' or 'LAT1'
# 'var4': 'VAR_NAME' or 'BAND2_VAL1'
# 'var5': 'LEVEL'
VARDICT={'var1':'YEAR',
         'var2':'MONTH',
         'var3':'BAND1_VAL1',
         'var4':'BAND2_VAL1',
         'var5':'LEVEL'}

######################################################
LOOPVAR1= ['X'] # dummy value if not needed
######################################################
#LOOPVAR1=[2019,]
#LOOPVAR1=range(1999,2018+1)

######################################################
LOOPVAR2= ['X'] # dummy value if not needed
######################################################
#LOOPVAR2=range(1,12+1) # Set MONTH ranges if outfile_frequency is less than 'year'
#LOOPVAR2=[3,]
#LOOPVAR2=['sg579m031oi01_zlev_h','sg534m031oi01_zlev_h','sg532m031oi01_zlev_h','sg620m031oi01_zlev_h','sg613m031oi01_zlev_h']
#LOOPVAR2=range(100) # NODE_NUMBER: Number of processors to farm out mean.py for Monte Carlo simulation of null distribution

######################################################
LOOPVAR3= ['X'] # dummy value if not needed
######################################################
#LOOPVAR3=['jan9815','feb9815','mar9815','apr9815','may9815','jun9815','jul9815','aug9815','sep9815','oct9815','nov9815','dec9815']
#LOOPVAR3=['ann'+str(xx).zfill(4) for xx in range(1979,2017+1)]
#LOOPVAR3=['ann9818','n2a9899-1718','m2o9818']
#LOOPVAR3=['-89.4629', '-88.7669', '-88.0669', '-87.3660', '-86.6648', '-85.9633', '-85.2618', '-84.5602', '-83.8586', '-83.1569', '-82.4553', '-81.7536', '-81.0519', '-80.3502', '-79.6485', '-78.9468', '-78.2450', '-77.5433', '-76.8416', '-76.1399', '-75.4381', '-74.7364', '-74.0346', '-73.3329', '-72.6312', '-71.9294', '-71.2277', '-70.5260', '-69.8242', '-69.1225', '-68.4207', '-67.7190', '-67.0172', '-66.3155', '-65.6137', '-64.9120', '-64.2102', '-63.5085', '-62.8067', '-62.1050', '-61.4033', '-60.7015', '-59.9998', '-59.2980', '-58.5963', '-57.8945', '-57.1928', '-56.4910', '-55.7893', '-55.0875', '-54.3858', '-53.6840', '-52.9823', '-52.2805', '-51.5788', '-50.8770', '-50.1753', '-49.4735', '-48.7718', '-48.0700', '-47.3683', '-46.6665', '-45.9647', '-45.2630', '-44.5612', '-43.8595', '-43.1577', '-42.4560', '-41.7542', '-41.0525', '-40.3507', '-39.6490', '-38.9472', '-38.2455', '-37.5437', '-36.8420', '-36.1402', '-35.4385', '-34.7367', '-34.0350', '-33.3332', '-32.6315', '-31.9297', '-31.2280', '-30.5262', '-29.8244', '-29.1227', '-28.4209', '-27.7192', '-27.0174', '-26.3157', '-25.6139', '-24.9122', '-24.2104', '-23.5087', '-22.8069', '-22.1052', '-21.4034', '-20.7017', '-19.9999', '-19.2982', '-18.5964', '-17.8947', '-17.1929', '-16.4911', '-15.7894', '-15.0876', '-14.3859', '-13.6841', '-12.9824', '-12.2806', '-11.5789', '-10.8771', '-10.1754', '-09.4736', '-08.7719', '-08.0701', '-07.3684', '-06.6666', '-05.9649', '-05.2631', '-04.5613', '-03.8596', '-03.1578', '-02.4561', '-01.7543', '-01.0526', '-00.3508', '00.3508', '01.0526', '01.7543', '02.4561', '03.1578', '03.8596', '04.5613', '05.2631', '05.9649', '06.6666', '07.3684', '08.0701', '08.7719', '09.4736', '10.1754', '10.8771', '11.5789', '12.2806', '12.9824', '13.6841', '14.3859', '15.0876', '15.7894', '16.4911', '17.1929', '17.8947', '18.5964', '19.2982', '19.9999', '20.7017', '21.4034', '22.1052', '22.8069', '23.5087', '24.2104', '24.9122', '25.6139', '26.3157', '27.0174', '27.7192', '28.4209', '29.1227', '29.8244', '30.5262', '31.2280', '31.9297', '32.6315', '33.3332', '34.0350', '34.7367', '35.4385', '36.1402', '36.8420', '37.5437', '38.2455', '38.9472', '39.6490', '40.3507', '41.0525', '41.7542', '42.4560', '43.1577', '43.8595', '44.5612', '45.2630', '45.9647', '46.6665', '47.3683', '48.0700', '48.7718', '49.4735', '50.1753', '50.8770', '51.5788', '52.2805', '52.9823', '53.6840', '54.3858', '55.0875', '55.7893', '56.4910', '57.1928', '57.8945', '58.5963', '59.2980', '59.9998', '60.7015', '61.4033', '62.1050', '62.8067', '63.5085', '64.2102', '64.9120', '65.6137', '66.3155', '67.0172', '67.7190', '68.4207', '69.1225', '69.8242', '70.5260', '71.2277', '71.9294', '72.6312', '73.3329', '74.0346', '74.7364', '75.4381', '76.1399', '76.8416', '77.5433', '78.2450', '78.9468', '79.6485', '80.3502', '81.0519', '81.7536', '82.4553', '83.1569', '83.8586', '84.5602', '85.2618', '85.9633', '86.6648', '87.3660', '88.0669', '88.7669', '89.4629']
#LOOPVAR3=['11.45','11.55']

######################################################
LOOPVAR4= ['X'] # dummy value if not needed
######################################################
#LOOPVAR4=['uwnd','vwnd','vrt','div','omega']
#LOOPVAR4=['uwnd','vwnd','vrt','div']
#LOOPVAR4=['dvrtdt','m_uwnd_dvrtdx','m_vwnd_dvrtdy','m_omega_dvrtdp','m_vrt_div','m_ff_div','m_beta_vwnd','m_domegadx_dvwnddp','domegady_duwnddp','source_dvrtdt','res_dvrtdt']
#LOOPVAR4=['dvrtdt','m_uwnd_dvrtdx','m_vwnd_dvrtdy','m_omega_dvrtdp','m_vrt_div','m_ff_div','m_beta_vwnd','m_domegadx_dvwnddp','domegady_duwnddp','source_dvrtdt','res_dvrtdt','vrt_horiz_adv','vrt_stretch','vrt_tilt']
#LOOPVAR4=['lat','lon','tsc','sa']
#LOOPVAR4=['121.45','121.55']

######################################################
#LOOPVAR5= ['X'] # dummy value if not needed
######################################################
#LOOPVAR5=[1000,975,950,925,900,875,850,825,800,775,750,700,650,600,550,500,450,400,350,300,250,225,200,175,150,125,100]
LOOPVAR5=[1000,925,850,700,600,500,400,300,250,200,150,100,70,50,30,20,10]

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
                    fout.write('#SBATCH -t 23:59:00\n')
                    fout.write('#SBATCH -p compute\n')
                    fout.write('#SBATCH --mem '+memory+'\n')
                    fout.write('#SBATCH --mail-type=NONE\n')
                    fout.write('#SBATCH --mail-user=e058@uea.ac.uk\n')
                    fout.write('#SBATCH --job-name='+FILENAME[:3]+'_'+str(index)+'\n')
                    fout.write('#SBATCH --nice='+str(NICE)+'\n')
                    if EXCLUSIVE:
                        fout.write('#SBATCH --exclusive\n')
                    fout.write('#SBATCH -o out_'+str(index)+'.txt\n')
                    fout.write('#SBATCH -e err_'+str(index)+'.txt\n')
                    fout.write('\n')
                    fout.write('python ./local_script.py\n')
                    fout.close()
            
                    # Change working directory to local directory
                    os.chdir(dir2)
                    
                    # Submit or run script
                    os.system('sbatch ./subfile')
                    
                    # Increment counter
                    index+=1

print('# Change back to initial working directory')
os.chdir(cwd)
