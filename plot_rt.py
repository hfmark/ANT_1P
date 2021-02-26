import numpy as np
import mods.PAT.utils as ut
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import ScientificColourMaps6 as SCM6
import os, sys

lon,lat,vph = np.loadtxt('rt_outputs/tomo_phase_3_8.1',unpack=True)
ext = [-77,-68,-54.5,-42]

plt.ion()
fig = plt.figure(figsize=(6.1,4.8))  # 9.2
ax1 = plt.subplot2grid((1,2),(0,0))  # 3
ax2 = plt.subplot2grid((1,2),(0,1))  # 3

# set up triangulation for grid
xi = np.linspace(min(lon),max(lon),20)
yi = np.linspace(min(lat),max(lat),30)
triang = tri.Triangulation(lon,lat)

# interpolate velocities
interp_vph = tri.LinearTriInterpolator(triang,vph)
Xi,Yi = np.meshgrid(xi,yi)
vph_grid = interp_vph(Xi,Yi)

# read and interpolate response amplitudes from resolution analysis
min_res = 0.02  # may need adjusting? but seems ok for now
lon,lat,rea = np.loadtxt('rt_outputs/tomo_phase_3_8.rea',usecols=(0,1,4),unpack=True)
interp_rea = tri.LinearTriInterpolator(triang,rea)
rea_grid = interp_rea(Xi,Yi)

lon,lat,dens = np.loadtxt('rt_outputs/tomo_phase_3_8.res',usecols=(0,1,2),unpack=True)
interp_dens = tri.LinearTriInterpolator(triang,dens)
dens_grid = interp_dens(Xi,Yi)

# mask velocities where no resolution
vph_grid.mask[rea_grid < min_res] = True

# plot masked array
con = ax1.imshow(vph_grid,extent=[min(xi),max(xi),min(yi),max(yi)],cmap=SCM6.roma,origin='lower')
ut.plot_basemap(ax1,bbox=ext)
cb = plt.colorbar(con,ax=ax1,fraction=0.06,orientation='horizontal')

con2 = ax2.imshow(dens_grid,extent=[min(xi),max(xi),min(yi),max(yi)],cmap=SCM6.hawaii,origin='lower')
ut.plot_basemap(ax2,bbox=ext)
cb2 = plt.colorbar(con2,ax=ax2,fraction=0.06,orientation='horizontal')
