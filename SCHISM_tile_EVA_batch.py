import os
import numpy as np
# from subprocess import call
import os

import glob
import re

sbatch_script_header = '''#!/bin/bash
#SBATCH --ntasks=4 #Number of CPU cores
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4 #Use other values here if using OpenMP instead of MPI
#SBATCH --mem=12GB
#SBATCH --nodes=1 #Number of nodes these cores should be shared between
#SBATCH -A OD-229130 #Must supply a currently active project code

# module load python
module load singularity
export SINGULARITYENV_PREPEND_PATH=/srv/conda/envs/notebook/bin:/srv/conda/condabin:/srv/conda/bin

'''

# scratch_path = '/scratch1/hoe01e/vankirap'
scratch_path = '/scratch3/hoe01e/vankirap'

# ncfiles = glob.glob(f'{scratch_path}/Hindcast_v2_hs_node_idxs_*.nc')

ncfiles = glob.glob(f'{scratch_path}/Hindcast_v3_hs_node_idxs_*.nc')
ncfiles.sort()

# for ncfile in [ncfiles[0]]:
for ncfile in ncfiles:
    node_idxs = re.sub('\.nc','', re.sub('Hindcast_v3_hs_node_idxs_','', ncfile.split('/')[-1]))
    outfile = f'{scratch_path}/Hindcast_v3_hs_EVA_node_idxs_{node_idxs}.nc'
    if os.path.isfile(outfile):
        continue
    else:
        node_idxs = re.sub('\.nc','', re.sub('Hindcast_v3_hs_node_idxs_','', ncfile.split('/')[-1]))
        sbatch_script = f'{scratch_path}/sch_EVA_{node_idxs}.bash'
        with open (sbatch_script, 'w') as rsh:
            rsh.write(sbatch_script_header)
            # rsh.writelines(f'python SCHISM_tile_EVA.py {ncfile}')
            rsh.writelines(f'singularity exec -B /home/hoe01e /datasets/work/oa-roamsurf/work/singularity/pangeo-airflow-code-20230615.sif python SCHISM_tile_EVA.py {ncfile}')
        print(f'submitting .... SCHISM_tile_EVA.py {ncfile}...')
        sbatch_call = f'sbatch {sbatch_script}'
        # call(sbatch_call, shell=False)
        os.system(sbatch_call)
