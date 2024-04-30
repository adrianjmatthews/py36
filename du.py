"""Disk usage of a given directory.

Uses unix du command. du itself does not return output in alphabetical order. This script does.
"""

import os


DIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data')
#DIR=os.path.join(os.path.sep,'gpfs','scratch','e058','data','imergv07amcw_sfc_30m')
#DIR=os.path.join(os.path.sep,'gpfs','afm','matthews','data')
#DIR=os.path.join(os.path.sep,'gpfs','afm','ocean_archive')

xx1=os.listdir(DIR)
xx1.sort()
for xx in xx1:
    #print(xx)
    command='du --block-size=1G --max-depth=0 '+os.path.join(DIR,xx)
    #print(command)
    os.system(command)

print('--------------------------')
command='du --block-size=1G --max-depth=0 '+os.path.join(DIR)
#print(command)
os.system(command)

