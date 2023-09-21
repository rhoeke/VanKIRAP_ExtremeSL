import os
import numpy as np
from subprocess import call

qsub_call = "qsub -l walltime=72:00:00 %s"
call(qsub_call % "test_job.sh", shell=True)


pbs_script_header = '''#!/bin/bash
#PBS -q normal
#PBS -l walltime=01:00:00
#PBS -l mem=64gb
#PBS -l ncpus=8
#PBS -l jobfs=100GB
#PBS -l wd
#PBS -j oe
#PBS -P v95
#PBS -l storage=gdata/hh5+gdata/v95+scratch/v95

module use /g/data/hh5/public/modules
module load conda/analysis3-22.10
'''

years = np.arange(1980, 2021)
# years = np.arange(1981, 1983)
output_path = '/scratch/v95/VanKIRAP/Hindcast_v3'
var2merge = {'elev':'elev','WWM_1':'hs'}
var_str = '_'.join(list(var2merge.values()))

for year in years:
    nc_merge_file = f'{output_path}/schout_{year}_{var_str}.nc'
    if os.path.exists(nc_merge_file):
        print(f'*{nc_merge_file} exists skipping!')
        continue
    else:
        pbs_script = f'sch_merge_{year}.sh'
        with open (pbs_script, 'w') as rsh:
            rsh.write(pbs_script_header)
            rsh.writelines(f'python hindcast_merge_year_files.py {year}')
        print(f'running .... hindcast_merge_year_files.py {year}...')
        qsub_call = f'qsub {pbs_script}'
        call(qsub_call, shell=True)
