import sys
# import os
import numpy as np
import xarray as xr
# import pandas as pd
import dask

''' This script merges (actually concatenates) selected variables of 
all SCHISM-WWM3 output files within a given (input) year

'''

year = sys.argv[1]

dask.config.set(scheduler='threads')

source_path = '/g/data/v95/VanKIRAP/UC/Hindcast_v3'

output_path = '/scratch/v95/VanKIRAP/Hindcast_v3'

years = np.arange(1980, 2021)

if not(int(year) in years):
    sys.exit(f'Input argument {year} is not a valid year!')

# file list in each year folder - this will fall over if the file is not there
# files schout_2.nc to schout_13.nc are Jan to Dec
# schout_1.nc is for spin-up, schout_14.nc is for next year's spin-up
nc_files= [f'schout_{i}.nc' for i in np.arange(2,14)]

# What variables? And rename them per key/value pair
var2merge = {'elev':'elev','WWM_1':'hs'}
# topologogy vars needed for plotting, etc.
topo_vars = ['SCHISM_hgrid_face_nodes','SCHISM_hgrid_node_x','SCHISM_hgrid_node_y',
        'depth','minimum_depth']

# Merge (concatenate) data variables (var2merge)
ds_list = []
for nc_file in nc_files:
    nc_path = f'{source_path}/{year}/outputs/{nc_file}'
    print(nc_path)
    ds = xr.load_dataset(nc_path)
    ds_subset=ds[list(var2merge.keys())]
    ds_subset = ds_subset.rename(var2merge)
    ds_list.append(ds_subset)

ds_merge = xr.concat(ds_list,dim='time')

# copy in topology variables from last dataset
ds_merge[topo_vars]=ds[topo_vars]


var_str = '_'.join(list(var2merge.values()))

nc_merge_file = f'{output_path}/schout_{year}_{var_str}.nc'
#nc_merge_file = f'{source_path}/schout_{year}_{var_str}.nc'
print(f'writing {nc_merge_file}')
ds_merge.to_netcdf(nc_merge_file)
