"""Loop over a variable and run scripts in series from an interactive session.

Useful for running short plotting scripts, as need connection to X
server which haven't figured out yet for batch jobs."""

import shutil

import datetime

import plot_tmpam_imerg_maps as plotscript

#LOOPVAR=range(27,31+1)

datec=datetime.datetime(2019,3,20)
date2=datetime.datetime(2019,4,14)
deltatime=datetime.timedelta(days=1)
LOOPVAR=[]
while datec<=date2:
    LOOPVAR.append(datec)
    datec=datec+deltatime

for loopvarc in LOOPVAR:
    print('######################################################')
    print(loopvarc)

    #shutil.copy(filec,tmpfile)
    #source_file=open(tmpfile)
    #exec(source_file.read())
    #del source_file
    #os.remove(tmpfile)

    #with open(filec) as source_file:
    #    exec(source_file.read())

    # Create instance of Plot object
    descriptor={}
    descriptor['MONTH']=loopvarc.month
    descriptor['DAY_START']=loopvarc.day
    aa=plotscript.Plot(**descriptor)
    print(aa)
    
    # Run plotting methods
    aa.run()
