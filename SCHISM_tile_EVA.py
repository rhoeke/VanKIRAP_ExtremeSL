import sys
import re
import numpy as np
import pandas as pd
import xarray as xr

from pyextremes import EVA
# from tqdm import tqdm


# from dask.distributed import Client, LocalCluster 
# cluster = LocalCluster() # Launches a scheduler and workers locally 
# client = Client(cluster)

ncfile = sys.argv[1]
node_idxs = re.sub('\.nc','', re.sub('Hindcast_v3_hs_node_idxs_','', ncfile.split('/')[-1]))


# scratch_path = '/scratch1/hoe01e/vankirap'
scratch_path = '/scratch3/hoe01e/vankirap'

da = xr.open_dataset(ncfile)
print(f'**** opened {ncfile}')
# da=da.isel(nSCHISM_hgrid_node=slice(0,20))
# print(f'**** reduced data {ncfile}')

eva_df=pd.DataFrame(index=da.nSCHISM_hgrid_node.values)
nnodes = eva_df.index.shape

return_period = [10, 50, 100]
eva_df['mean'] = da.hs.chunk(dict(time=-1)).mean(dim='time').values
eva_df['max'] = da.hs.chunk(dict(time=-1)).max(dim='time').values
eva_df['p99'] = da.hs.chunk(dict(time=-1)).quantile(0.99, dim='time').values
dist_params = ['c', 'scale']
nan_col = np.ones(nnodes)*np.nan

for col in dist_params:
    eva_df[col]=nan_col
for rt in return_period:
    eva_df[f'{rt}']=nan_col
    eva_df[f'{rt} lower ci']=nan_col
    eva_df[f'{rt} upper ci']=nan_col

# for idx in eva_df.index:
# for idx in tqdm(range(nnodes[0]), desc='Iterations completed'):
for idx in np.arange(0,nnodes[0]): 
    if eva_df['p99'].loc[idx]>0.25:
        try:
            data = da.hs.isel(nSCHISM_hgrid_node=idx).to_series()
            model = EVA(data)
            thres = eva_df['p99'].loc[idx]
            model.get_extremes(method="POT", threshold=thres, r="24H")
            model.fit_model()
            summary = model.get_summary(return_period=return_period, alpha=0.95)
            for key, val in model.distribution.mle_parameters.items():
                eva_df.loc[idx,key] = val
            for rt in return_period:
                eva_df.loc[idx,f'{rt}'] = summary.loc[rt,'return value']
                eva_df.loc[idx,f'{rt} lower ci'] = summary.loc[rt,'lower ci']
                eva_df.loc[idx,f'{rt} upper ci'] = summary.loc[rt,'upper ci']
        except:
            print('*** something is wrong with idx={idx} in node_idxs')
print('**** Completed EVA calc, converting to dataset and saving!')
ds = xr.Dataset.from_dataframe(eva_df)
ds.to_netcdf(f'{scratch_path}/Hindcast_v3_hs_EVA_node_idxs_{node_idxs}.nc')
