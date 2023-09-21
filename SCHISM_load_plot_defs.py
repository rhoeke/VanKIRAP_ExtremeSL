import os
import re
import numpy as np
import xarray as xr
import pandas as pd 

from pyproj import Proj, transform

import matplotlib.pyplot as plt
# from mpl_toolkits.basemap import Basemap
from matplotlib.tri import Triangulation

# use cartopy rather than defunct basemap
import cartopy.crs as ccrs
import cartopy.feature as cfeature
#import cartopy.mpl.ticker as cticker
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature

from scipy.spatial import KDTree
from scipy.interpolate import griddata
'''
These functions have been copied from:
https://bitbucket.csiro.au/projects/CMEXTREMES/repos/ppbcha/browse/scripts/SCHISM/SCHISM_output_plotting_funcs.py

'''

def schism_meshtri(schout_ds):
    ''' 
    creates matplotlib trimesh object for plotting
    args: sc_ds: xarray dataset representation of SCHISM output file
    '''
    # get elements
    # For some reason, the SCHISM elements variable ("face_nodes") uses 1-based indexing, 
    # and has a fourth column of NaNs (maybe for 3D meshes?). Also the indices are stored as floats
    elems = int(schout_ds.SCHISM_hgrid_face_nodes[:,:-1]-1)
    # get lat/lon coordinates of nodes - weird it appears x,y are switched 
    lons = schout_ds.SCHISM_hgrid_node_y.values
    lats = schout_ds.SCHISM_hgrid_node_x.values
    # create trimesh object
    meshtri = Triangulation(lats, lons, triangles=elems)
    return meshtri

def schism_load(schfile):
    ''' 
    Load schism output file and parse the element/node values to create a (matplotlib) trimesh object for plotting
    
    Returns xarray dataset and trimesh object
    We should modify this to load multiple files ... probably need assistance from DASK
    '''
    schout_ds = xr.open_dataset(schfile)
    # get elements
    # For some reason, the SCHISM elements variable ("face_nodes") uses 1-based indexing, 
    # and has a fourth column of NaNs (maybe for 3D meshes?). Also the indices are stored as floats
    elems = int(schout_ds.SCHISM_hgrid_face_nodes[:,:-1]-1)
    # get lat/lon coordinates of nodes - weird it appears x,y are switched 
    lons = schout_ds.SCHISM_hgrid_node_y.values
    lats = schout_ds.SCHISM_hgrid_node_x.values
    # create trimesh object
    meshtri = Triangulation(lats, lons, triangles=elems)
    return schout_ds, meshtri


#define functions
#Cartesian convention
def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return x, y
#Nautical convention
def pol2cart2(rho, deg):
    x, y = pol2cart(rho, deg/180.*np.pi)
    return y, x

def near2d_unstruct(x, y, x0, y0):
    """
    Find the indexes of the mesh node that is
    nearest a chosen (x0, y0).
    Usage: line, col = near2d(x, y, x0, y0)
    
    **** should replace this using KDTree!
    
    """   
    
    
    dx = np.abs(x - x0); dx = dx / dx.max()
    dy = np.abs(y - y0); dy = dy / dy.max()
    dn = dx + dy    

    fn = np.where(dn==dn.min())[0][0]
    return int(fn)

def find_nearest_wet_point(schout,lon,lat,depth_thresh=0):
    '''Do an expanding search for nearest kn points untill the nearest point with depths>threshold is found
    '''
    # Grid Depths
    grid_depths=schout.depth.values
    # schism depths are by definition positive values (land is negative) so make -z down:
    grid_depths=grid_depths*-1
    
    # Project lons/lat into UTM coords 
    # - note this assumes that lon/lat are in WGS84 and we are projecting into UTM Zone 55s/WGS84 - this was for Vic
    # pUTM = Proj('epsg:32755')
    # - note this assumes that lon/lat are in WGS84 and we are projecting into UTM Zone 59s/WGS84 - this is for Vanuatu
    pUTM = Proj('epsg:32759')
    x,y = pUTM(schout.SCHISM_hgrid_node_x.values,schout.SCHISM_hgrid_node_y.values)
    x0,y0=pUTM(lon,lat)
    # Build KDtree
    tree=KDTree(list(zip(x,y)))
    input_point=list((float(x0),float(y0)))
    
    # first go - start with nearest single node
    dist,grid_idx=tree.query(input_point, k=1)
    depths=grid_depths[grid_idx]
    # successive increase by  5 nearest neighboors
    kn = 5
    while depths.min()>=depth_thresh:
        dist,grid_idx=tree.query(input_point, k=kn)
        depths=grid_depths[grid_idx]
        # increase nearest neighbor search by 5 points each loop
        kn=kn+5
#         print(kn)
        if dist.max()>3000: 
            print('*** warning search radius >3 km (%3.1f m, %d points), exiting search!'%(dist.max(),kn))
            break
    # This only needs to run if the nearest single node depth is not less than the threshold
    if grid_idx.size>1:
        idx2 = np.where(depths<=depth_thresh)
        grid_idx = grid_idx[idx2]
        depths = depths[idx2]
        depth_idx=np.where(depths==depths.max())[0]
        if len(depth_idx)>1:
            print(f'*** Warning, more than one point found with equally deep depths (z = {depths.max()}), defauling to first point')
            depth_idx=depth_idx[0]
        grid_idx=grid_idx[depth_idx]
    return int(grid_idx),kn



def schism_calculate_vectors(ax, schout, vtype='waves', dX='auto', mask=True):
    pUTM55 = Proj('epsg:32755')
    # pWGS84 = Proj('epsg:4326')
    if vtype=='waves':
        idx=(schout.WWM_1>0.05) & (schout.elev-schout.depth<0.1)
        dp=schout.WWM_18[idx]
        # hs=schout.WWM_1[idx]
        hs=np.ones(dp.shape)
        [u,v] = pol2cart2(hs, np.mod(dp+180, 360))
    elif vtype=='elev'or re.search('curr',vtype):
        idx=np.sqrt(schout.dahv[:,0]**2+schout.dahv[:,1]**2)>0.01
        u = schout.dahv[idx,0] 
        v = schout.dahv[idx,1] #dahv has u and v components, so use index of 1 for v and index of 0 for u
    else: 
        print('*** Warning input vecter data not understood')
    x,y = pUTM55(schout.SCHISM_hgrid_node_x.values[idx],schout.SCHISM_hgrid_node_y.values[idx])
    xlim,ylim=pUTM55(ax.get_xlim(),ax.get_ylim())
    # might have to play with this - we assume a total of 100 arrows a side will be pleasing
    if dX=='auto':
        n=30
        dX = np.ceil((ylim[1]-ylim[0])/n)
    xi = np.arange(xlim[0],xlim[1]+dX,dX)
    yi = np.arange(ylim[0],ylim[1]+dX,dX)
    XI,YI = np.meshgrid(xi,yi)

    UI = griddata((x,y),u,(XI,YI),method='linear')
    VI = griddata((x,y),v,(XI,YI),method='linear')
    # Create a mask so that place with very little data are removed
    if mask:
        xbins = np.arange(xlim[0],xlim[1]+2*dX,dX)
        ybins = np.arange(ylim[0],ylim[1]+2*dX,dX)
        densityH,_,_ = np.histogram2d(x, y, bins=[xbins,ybins])
        densityH=densityH.T
        # might want to adjust this...
        idx=densityH<1
        UI[idx]=np.NaN
        VI[idx]=np.NaN
    LonI,LatI = pUTM55(XI,YI, inverse=True)

    return LonI,LatI,UI,VI

def schism_plot(schout, meshtri,varname,varscale=[],bbox=[],time=-1,mask=True,
                vectors=False, plotmesh=False,project=False,contours=[10,30,50], 
                pscale=20,cmap=plt.cm.jet,cblabel = ''):
    ''' 
    plot output variable in xarray dataset (schout) using mesh information meshtri.
    input:
        schout: xarray dataset returned by def schism_load
        meshtri: (matplotlib) trimesh object returned by def schism_load
        varname: name of data variable in schout
        OPTIONAL:
            varscale: min/max plotting colour values (defalts to min/max)
            bbox: bounding box for plotting [minlon,minlat,maxlon,maxlat] (defalts to data bounds)
            time: time to plot (if a variable dimension), can be int (index) or datetime64-like object
            plotmesh: plot the grid mesh (elemment boundaries)
            project: use a map projection (cartopy, so that e.g. gis data can be added - this is slower
            mask: mask out parts of the SCHISM output based on a minumum depth threshold
    Returns xarray dataset and 
    We should modify this to load multiple files ... probably need assistance from DASK
    '''
    if 'time' in list(schout.dims):
        if type(time)==int : # input ts is index
            schout=schout.isel(time=time)
        else: # assume is datetime64-like object
            schout=schout.sel(time=time)
    # By default, plot depth contours 
    # ... I like depth to be z (negative)
    z=schout.depth*-1
    if np.mean(contours)>0:
        contours = np.flip(-1*np.asarray(contours))
    else:
        contours = np.asarray(contours)
    if varname=='depth' or varname=='z':
        var = z
    else:
        var=schout[varname]

    if len(varscale)==0: 
        vmin=var.min()
        vmax=var.max()
    else:
        vmin,vmax=varscale
    if project:
        x,y = meshtri.x,meshtri.y
        fig, ax = plt.subplots(1, 1, figsize=(pscale,pscale),
            subplot_kw={'projection': ccrs.PlateCarree()})
        if len(bbox)==4:
            ax.set_extent([bbox[0], bbox[2], bbox[1], bbox[3]], ccrs.PlateCarree())
        else:
            bbox=[x.min(),y.min(),x.max(),y.max()]
    else:
        fig, ax = plt.subplots(1, 1, figsize=(30,30))
    if plotmesh:
        ax.triplot(meshtri, color='k', alpha=0.3)
    ### MASKING ** doesnt work with tripcolor, must use tricontouf ###############################
    # mask all places in the grid where depth is greater than elev (i.e. are "dry") by threshold below
    if mask:
        # temp=var.values
        # threshold of + 0.05 seems pretty good   *** But do we want to use the minimum depth
        # defined in the SCHISM input (H0) and in output schout.minimum_depth
        # bad_idx= schout.elev.values+schout.depth.values<0.05
        bad_idx= schout.elev.values+schout.depth.values < schout.minimum_depth.values
        # new way
        mask = np.all(np.where(bad_idx[meshtri.triangles], True, False), axis=1)
        meshtri.set_mask(mask)
        #old way:
        # temp[bad_idx]=np.NaN
        # temp = np.ma.masked_invalid(temp)
        # levels = np.linspace(vmin, vmax, 50)
        # temp = temp.filled(fill_value=-999)
        # var.values=temp
        extend='neither'
        if (var.min()<vmin) & (var.max()>vmax):
            extend='both'
        elif var.min()<vmin:
            extend='min'
        elif var.max()>vmax:
            extend='max'

        cax = ax.tricontourf(meshtri, var, cmap=cmap,levels=np.linspace(vmin, vmax, 50), extend=extend)
    #no mask#############################################################
    else:
        cax = ax.tripcolor(meshtri, var, cmap=cmap, vmin=vmin, vmax=vmax)
        # quiver variables if asked
    if vectors:
        if re.search('WWM',varname):
            vtype='waves'
        else:
            vtype='currents'
        LonI,LatI,UI,VI=schism_calculate_vectors(ax, schout, vtype=vtype)
        ax.quiver(LonI,LatI,UI,VI, color='k')

    con = ax.tricontour(meshtri, z, contours, colors='k')
    # ax.clabel(con, con.levels, inline=True, fmt='%i', fontsize=12)
    if not(project):
        ax.set_aspect('equal')
    else:
        # draw lat/lon grid lines every n degrees.
        #  n = (max(lat)-min(lat))/8
        n = (bbox[2]-bbox[0])/5
        for fac in [1,10,100]:
            nr = np.round(n*fac)/fac
            if nr>0:
                n=nr
                xticks = np.arange(np.round(bbox[0]*fac)/fac,np.round(bbox[2]*fac)/fac+n,n)
                yticks = np.arange(np.round(bbox[1]*fac)/fac,np.round(bbox[3]*fac)/fac+n,n)
                break
#        ax.set_xticks(xticks, crs=ccrs.PlateCarree()
        ax.set_xticks(xticks)
#         ax.set_xticklabels(np.arange(np.round(min(x)),np.round(max(x)),n))
#        ax.set_yticks(yticks, crs=ccrs.PlateCarree()
        ax.set_yticks(yticks)
#         ax.set_yticklabels(np.arange(np.round(min(y)),np.round(max(y)),n))
        ax.yaxis.tick_left()
        #lon_formatter = cticker.LongitudeFormatter()
        #lat_formatter = cticker.LatitudeFormatter()
        # New versions of marplotlib throw warnings on this - does it matter
	# ax.xaxis.set_major_formatter(lon_formatter)
        # ax.yaxis.set_major_formatter(lat_formatter)
	#ax.set_xticks(lon_formatter)
	#ax.set_yticks(lat_formatter)
        ax.grid(linewidth=1, color='black', alpha=0.5, linestyle='--')
        ax.add_feature(cfeature.BORDERS, linewidth=2)
    if len(bbox)==4:
        ax.set_ylim(bbox[1], bbox[3])
        ax.set_xlim(bbox[0], bbox[2])

    cbar = plt.colorbar(mappable=cax,shrink=0.5)
    cbar.set_ticks(np.round(np.linspace(vmin,vmax,5)*100)/100)
    if cblabel=='':
        cbar.set_label(varname)
    else:
        cbar.set_label(cblabel)
    return fig, ax