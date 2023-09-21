# This collection of Python notebooks and scripts/functions is to perform extreme value analysis (EVA) on the VAnKIRAP SCHISM (deterministic) Hindcast
Most things are assumed to be run on CSIRO HPC/Bowen cloud resources so some paths are hard-coded; however some initial processing of SCHISM output on Gadi/VDI (NCI) is required (see next section).

## NCI SCHISM national hindcast (CFSR/CAWCR/ORAS5) output data preparation:
The "raw" model files are here (NCI):
`/g/data/v95/VanKIRAP/UC/Hindcast_v3/`

The individual parallel processed output have (hopefully) been "merged" into monthly outputs.  There are then two simple python script (in this repo) needed to futher concanate a subset of variables in preparation for transfer to CSIRO infrastruction:
1. `hindcast_merge_year_files.py`: selects and renames variables (currently {'elev':'elev','WWM_1':'hs'}) and merges them in single yearly files, saving them to `/scratch/v95/VanKIRAP/`
2. `submit_hindcast_merge_jobs.py`: a wrapper for the above script, running all years, each as a PBS job
The output of these two looks like:
schout_1981_elev_hs.nc  
schout_1982_elev_hs.nc
...

## SCHISM data EVA analysis on CSIRO infrastructure:

This presumes that the yearly, subset merged files processed per the above section have been moved/compied to a location visible to (e.g.) Petrichor:
'ls /datasets/work/oa-vankirap/work/schism/Hindcast_v3/
schout_1981_elev_hs.nc  
schout_1982_elev_hs.nc
...'

The remainder of this section can be divided into two sub-sections, *EVA of coastal points* and *EVA of the entire native SCHISM mesh*

### EVA of coastal points

*Extract coastal points workflow:
1. `Coastal_points_gen.ipynb`:  reads a coastline shapefile and creates regularly spaced coastal point locations (saved as separate shapefile) for EVA
2. `SCHISM_check_problem_areas.ipynb`: checking/plotting areas of insuffiecient resolution and/or wetting-drying problems in the mesh output
3. `SCHISM_extract_coastal_point_data.ipynb`: reads coastpoints shapefile and finds nearest logical "wet-point" in SCHISM Mesh, exacts relavent time-series variable (e.g. 'elev') from SCHISM hindcast dataset, saving as new nc file
4. `SCHISM_coastal_points_EVA.ipynb`: perform EVA analysis (including calculation of 10,50 and 100-year ARIs) on coastal point time series from previous step and save output files (e.g. to provide to VanKIRAP portal).  Uses pyextremes package, but could use something else ... 
5. `SCHISM_coastal_points_EVA_add_SLR.ipynb`: lineararly add in (Zhang et al 2021) sea level rise scenarios to 10,50 and 100-year ARIs ..


**Note 1** Step 4 and 5 outputs the coastal points EVA (ARIs) as csv files - we should probably fix this up to output netCDF with associated metadata ...


### EVA of the native SCHISM mesh

1. `SCHISM_EVA_waves_tile.ipynb`:
2. sbatch/SLURM proccessing:
	- `SCHISM_tile_EVA.py`
    - `SCHISM_tile_EVA_batch.py`
3. `SCHISM_EVA_waves_tile.ipynb`:

