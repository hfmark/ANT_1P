import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import mods.PAT.anttomo as ant
import mods.PAT.files as vr
import cartopy.crs as ccrs
import ScientificColourMaps6 as SCM6
import pickle
import os, sys

####
# read maps from EQ and ANT tomo, compare at one period. ISOTROPIC.
####

# read in EQ tomography outputs
#aphv = sio.loadmat(os.path.expanduser('~/Patagonia/EQ_tomo/helmholtz_stack_ENAPnoamp_noMG.mat'),variable_names='avgphv')
#aphv = sio.loadmat(os.path.expanduser('~/Patagonia/EQ_tomo/helmholtz_stack_ENAP_noMG_noRMG.mat'),variable_names='avgphv')
aphv = sio.loadmat(os.path.expanduser('~/Patagonia/EQ_tomo/helmholtz_stack_all.mat'),variable_names='avgphv')
aphv = aphv['avgphv']

periods = [20,26,32,40]
ip = 0

# pick a period to compare
if len(aphv[0]) == 8:  # just the test periods
    einds = [0,1,2]
elif len(aphv[0]) == 24:  # all the periods
    einds = [0,3,6,10]

period = periods[ip]
eind = einds[ip]
vlm = vr.vlims[period]

#aphv = ant.EQVelocityMap(aphv[:,eind],eik=True)
aphv = ant.EQVelocityMap(aphv[:,eind],eik=False)

# read in some station info
slon, slat = np.loadtxt(vr.sta_list,usecols=(1,2),unpack=True)

# start plotting
plt.ion()
fig = plt.figure(figsize=(9.2,4.8))
ax1 = plt.subplot2grid((1,3),(0,0),projection=ccrs.PlateCarree())
ax2 = plt.subplot2grid((1,3),(0,1),projection=ccrs.PlateCarree())
ax3 = plt.subplot2grid((1,3),(0,2),projection=ccrs.PlateCarree())
xtk = np.arange(-75,-68,3)
ytk = np.arange(-55,-43,3)

im1 = ax1.imshow(aphv.v.T,\
		extent=aphv.grid.bbox(),transform=ccrs.PlateCarree(),\
		origin='lower',cmap=SCM6.roma)  # EQ tomo
im1.set_clim(vlm)

# read raytomo outputs
lon,lat,vph = np.loadtxt('rt_outputs/tomo_iso_phase_3_%i.1' % period, usecols=(0,1,2),unpack=True)

# set up triangulation for grid
xi = np.arange(min(lon),max(lon),0.3)
yi = np.arange(min(lat),max(lat),0.3)

# interpolate velocities
triang = tri.Triangulation(lon,lat)
interp_vph = tri.LinearTriInterpolator(triang,vph)
Xi,Yi = np.meshgrid(xi,yi)
vph_grid = interp_vph(Xi,Yi)

# read and interpolate response amplitudes from resolution analysis
min_res = 0.05  # may need adjusting? but seems ok for now
lon,lat,rea = np.loadtxt('rt_outputs/tomo_iso_phase_3_%i.rea' % (period),usecols=(0,1,4),unpack=True)
triang = tri.Triangulation(lon,lat)
interp_rea = tri.LinearTriInterpolator(triang,rea)
rea_grid = interp_rea(Xi,Yi)

# read and interpolate orbit density
min_dens = 10 # may need adjusting?
lon,lat,dens = np.loadtxt('rt_outputs/tomo_iso_phase_3_%i.res' % (period),usecols=(0,1,2),unpack=True)
triang = tri.Triangulation(lon,lat)
interp_dens = tri.LinearTriInterpolator(triang,dens)
dens_grid = interp_dens(Xi,Yi)

# mask velocities where no resolution and/or not enough paths
vph_grid.mask[rea_grid < min_res] = True
vph_grid.mask[dens_grid < min_dens] = True

# pad grid to make it match the eq tomo [NOTE this assumes certain shapes]
to_add = vph_grid.shape[1]
vph_grid = np.ma.vstack((vph_grid,np.ma.MaskedArray(np.zeros(to_add),mask=np.ones(to_add))))

# plot isotropic raytomo
im2 = ax2.imshow(vph_grid,\
                extent=[min(xi),max(xi),min(yi),max(yi)],transform=ccrs.PlateCarree(),\
		origin='lower',cmap=SCM6.roma)
im2.set_clim(vlm)

# difference the grids
diff = np.zeros(aphv.v.shape)
# loop grid of EQ map (aphv); if nan, add nan; if not nan, interpolate ANT and subtract
lons,lats = aphv.grid.xy_nodes()
for i in range(len(lons)):
    lon = lons[i]; lat = lats[i]
    ix,iy = aphv.grid.ix_iy(i)
    if np.isfinite(aphv.v[ix,iy]) and np.isfinite(vph_grid.T[ix,iy]): # NOTE: assumes same gridding
        avel = vph_grid.T[ix,iy]
        diff[ix,iy] = avel - aphv.v[ix,iy]
    else:
        diff[ix,iy] = np.nan

im3 = ax3.imshow(diff.T,\
		extent=aphv.grid.bbox(),transform=ccrs.PlateCarree(),\
		origin='lower',cmap=SCM6.cork,vmin=-0.2,vmax=0.2)
ax3.contour(diff.T,levels=(-0.1,0.1),colors=('r','m'),extent=aphv.grid.bbox(),origin='lower')

# general plot things
#ant.basemap(ax=ax1,bbox=aphv.grid.bbox())
#ant.basemap(ax=ax2,bbox=aphv.grid.bbox())
#ant.basemap(ax=ax3,bbox=aphv.grid.bbox())
for ax in [ax1,ax2,ax3]:
    ax.coastlines(resolution='50m')
    ax.gridlines(crs=ccrs.PlateCarree(),xlocs=xtk,ylocs=ytk)
    ax.set_xticks(xtk,minor=False,crs=ccrs.PlateCarree())
    ax.set_yticks(ytk,minor=False,crs=ccrs.PlateCarree())

plt.suptitle('%.1f s' % period)
ax3.yaxis.tick_right()
ax3.yaxis.set_label_position('right')
for tl in ax2.get_yticklabels():
    tl.set_visible(False)
ax2.set_ylabel('')

cb1 = plt.colorbar(im1,ax=ax1,fraction=0.04,orientation='horizontal')                       
cb2 = plt.colorbar(im2,ax=ax2,fraction=0.04,orientation='horizontal')                       
cb3 = plt.colorbar(im3,ax=ax3,fraction=0.04,orientation='horizontal')


sys.exit()






ax1.scatter(slon,slat,s=20,marker='^',color='k')
ax2.scatter(slon,slat,s=20,marker='^',color='k')
ax3.scatter(slon,slat,s=20,marker='^',color='k')

                      

