import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import mods.PAT.utils as ut
import mods.PAT.files as vr
import ScientificColourMaps6 as SCM6
import pickle
import os, sys

####
# read maps from anttomo and RayTomo and compare at one period
####

period = 32.

# set up figure
plt.ion()
fig = plt.figure(figsize=(12.3,4.8))
ax1 = plt.subplot2grid((1,4),(0,0))
ax2 = plt.subplot2grid((1,4),(0,1))
ax3 = plt.subplot2grid((1,4),(0,2))
ax4 = plt.subplot2grid((1,4),(0,3))

# read anttomo outputs
f = open('output/4-pass-tomography_phase.pickle','rb')
vmaps = pickle.load(f)
f.close()

try:
    vmaps[period]
except KeyError:
    print('missing period in anttomo output')
    sys.exit()

# get anttomo velocities
vant = vmaps[period].grid.to_2D_array(vmaps[period].v0 / (1 + vmaps[period].mopt))
amap = vmaps[period]

# put nans in any cell with fewer than 10 paths
dens = amap.grid.to_2D_array(amap.density)
for i in range(dens.shape[0]):
    for j in range(dens.shape[1]):
        if dens[i,j] < 10:
            vant[i,j] = np.nan

# plot anttomo
im1 = ax1.imshow(vant.T,\
		extent=amap.grid.bbox(),
		origin='lower',cmap=SCM6.roma)
im1.set_clim(vr.vlims[period])

# read isotropic RayTomo outputs
lon,lat,vph = np.loadtxt('rt_outputs/tomo_iso_phase_3_%i.1' % (period),usecols=(0,1,2),unpack=True)

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

# plot isotropic RayTomo
im2 = ax2.imshow(vph_grid,\
                extent=[min(xi),max(xi),min(yi),max(yi)],\
		origin='lower',cmap=SCM6.roma)
im2.set_clim(vr.vlims[period])


# read anisotropic RayTomo outputs
lon,lat,vph = np.loadtxt('rt_outputs/tomo_ani_phase_3_%i.1' % (period),usecols=(0,1,2),unpack=True)

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
lon,lat,rea = np.loadtxt('rt_outputs/tomo_ani_phase_3_%i.rea' % (period),usecols=(0,1,4),unpack=True)
triang = tri.Triangulation(lon,lat)
interp_rea = tri.LinearTriInterpolator(triang,rea)
rea_grid = interp_rea(Xi,Yi)

# read and interpolate orbit density
min_dens = 10 # may need adjusting?
lon,lat,dens = np.loadtxt('rt_outputs/tomo_ani_phase_3_%i.res' % (period),usecols=(0,1,2),unpack=True)
triang = tri.Triangulation(lon,lat)
interp_dens = tri.LinearTriInterpolator(triang,dens)
dens_grid = interp_dens(Xi,Yi)

# mask velocities where no resolution and/or not enough paths
vph_grid.mask[rea_grid < min_res] = True
vph_grid.mask[dens_grid < min_dens] = True

# plot anisotropic RayTomo
im3 = ax3.imshow(vph_grid,\
                extent=[min(xi),max(xi),min(yi),max(yi)],\
		origin='lower',cmap=SCM6.roma)
im3.set_clim(vr.vlims[period])


# plot anisotropic RayTomo with arrows
im4 = ax4.imshow(vph_grid,\
                extent=[min(xi),max(xi),min(yi),max(yi)],\
		origin='lower',cmap=SCM6.roma)
im4.set_clim(vr.vlims[period])
lon,lat,a_amp,a_psi = np.loadtxt('rt_outputs/tomo_ani_phase_3_%i.1' % (period),\
            usecols=(0,1,5,6),unpack=True)
ax4.quiver(lon,lat,a_amp*np.cos(np.degrees(a_psi)),a_amp*np.sin(np.degrees(a_psi)))

################################################3
# fancy up the plots
cb1 = plt.colorbar(im1,ax=ax1,fraction=0.04,orientation='horizontal')
cb2 = plt.colorbar(im2,ax=ax2,fraction=0.04,orientation='horizontal')
cb3 = plt.colorbar(im3,ax=ax3,fraction=0.04,orientation='horizontal')
cb4 = plt.colorbar(im4,ax=ax4,fraction=0.04,orientation='horizontal')

ut.plot_basemap(ax1,bbox=amap.grid.bbox())
ut.plot_basemap(ax2,bbox=amap.grid.bbox())
ut.plot_basemap(ax3,bbox=amap.grid.bbox())
ut.plot_basemap(ax4,bbox=amap.grid.bbox())

ax1.set_title('python code')
ax2.set_title('rt isotropic')
ax3.set_title('rt anisotropic')

plt.suptitle('%i sec' % period)
