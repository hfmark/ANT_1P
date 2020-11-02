import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import mods.PAT.anttomo as ant
import mods.PAT.files as vr
import pickle
import os, sys

####
# read maps from EQ and ANT tomo, compare at one period
####

# read in EQ tomography outputs
aphv = sio.loadmat(os.path.expanduser('~/Patagonia/EQ_tomo/save_noENAP/helmholtz_stack_all.mat'),variable_names='avgphv')
aphv = aphv['avgphv']

# pick a period to compare
if len(aphv[0]) == 8:  # just the test periods
    #eind = 0; aind = 20.0; vmin=3.22; vmax=3.8  # eind for test period file
    eind = 1; aind = 25.0; vmin = 3.429; vmax = 3.98
    #eind = 2; aind = 32.0; vmin = 3.43; vmax = 4.07
elif len(aphv[0]) == 24:  # all the periods
    #eind = 0; aind = 20.0; vmin=3.22; vmax = 3.8
    eind = 3; aind = 26.0; vmin = 3.429; vmax = 3.98
    #eind = 6; aind = 32.0; vmin = 3.43; vmax = 4.07


#aphv = ant.EQVelocityMap(aphv[:,eind],eik=True)
aphv = ant.EQVelocityMap(aphv[:,eind],eik=False)

# read in ANT outputs
#f = open('output/3-pass-tomography_phase.pickle','rb')
f = open('output/4-pass-tomography_phase.pickle','rb')
vmaps = pickle.load(f)
f.close()
try:
    vmaps[aind]
except IndexError:
    aind = 25.0  # 25/26 potential mismatch if EQ is all and ANT is test
vant = vmaps[aind].grid.to_2D_array(vmaps[aind].v0 / (1 + vmaps[aind].mopt))
amap = vmaps[aind]

# put nans in any cell with fewer than 10 paths
dens = amap.grid.to_2D_array(amap.density)
for i in range(dens.shape[0]):
    for j in range(dens.shape[1]):
        if dens[i,j] < 10:
            vant[i,j] = np.nan

# calculate difference at EQ grid points
diff = np.zeros(aphv.v.shape)
# loop grid of EQ map (aphv); if nan, add nan; if not nan, interpolate ANT and subtract
lons,lats = aphv.grid.xy_nodes()
for i in range(len(lons)):
    lon = lons[i]; lat = lats[i]
    ix,iy = aphv.grid.ix_iy(i)
    if np.isfinite(aphv.v[ix,iy]) and np.isfinite(vant[ix,iy]): # NOTE: assumes same gridding
        avel = amap.interp_velocity(lon,lat)
        diff[ix,iy] = avel - aphv.v[ix,iy]
    else:
        diff[ix,iy] = np.nan

# set vmin/vmax
#m1 = np.nanmin(aphv.v)
#m2 = np.nanmin(vant)
#vmin = min(m1,m2)
#m1 = np.nanmax(aphv.v)
#m2 = np.nanmax(vant)
#vmax = max(m1,m2)


# read in some station info
slon, slat = np.loadtxt(vr.sta_list,usecols=(1,2),unpack=True)

plt.ion()
fig = plt.figure(figsize=(9.2,4.8))
ax1 = plt.subplot2grid((1,3),(0,0))
ax2 = plt.subplot2grid((1,3),(0,1))
ax3 = plt.subplot2grid((1,3),(0,2))
im1 = ax1.imshow(aphv.v.T,\
		extent=amap.grid.bbox(),\
		origin='bottom',cmap='seismic_r',vmin=vmin,vmax=vmax)

im2 = ax2.imshow(vant.T,\
		extent=amap.grid.bbox(),
		origin='bottom',cmap='seismic_r',vmin=vmin,vmax=vmax) #,\
	       # interpolation='bicubic')

im3 = ax3.imshow(diff.T,\
		extent=amap.grid.bbox(),
		origin='bottom',cmap='PuOr')

#r = amap.grid.to_2D_array(amap.Rradius)
#m1 = ax2.imshow(r.T,\
#	        origin='bottom', extent=amap.grid.bbox(),\
#		interpolation='bicubic',cmap=ant.CMAP_MASK)  

ant.basemap(ax=ax1,bbox=amap.grid.bbox())
ant.basemap(ax=ax2,bbox=amap.grid.bbox())
ant.basemap(ax=ax3,bbox=amap.grid.bbox())

plt.suptitle('%.1f s' % aind)
ax3.yaxis.tick_right()
ax3.yaxis.set_label_position('right')
for tl in ax2.get_yticklabels():
    tl.set_visible(False)
ax2.set_ylabel('')

ax1.scatter(slon,slat,s=20,marker='^',color='k')
ax2.scatter(slon,slat,s=20,marker='^',color='k')
ax3.scatter(slon,slat,s=20,marker='^',color='k')

cb1 = plt.colorbar(im1,ax=ax1,fraction=0.04,orientation='horizontal')                       
cb2 = plt.colorbar(im2,ax=ax2,fraction=0.04,orientation='horizontal')                       
cb3 = plt.colorbar(im3,ax=ax3,fraction=0.04,orientation='horizontal')                       

