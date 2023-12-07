#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 16:26:26 2019

@author: fiercenator
"""

import numpy as np
import os
import netCDF4
import matplotlib.pyplot as plt
# from mpl_toolkits.basemap import Basemap
from basemap import Basemap
from matplotlib.table import table
from matplotlib import rc

read_dir = '/Users/fier887/Downloads/02_aging_timescale/gcm_data/CAM_aging/'
toolbox = '../../../process_partmc/'

# toolbox = '../../code/process_partmc/'
cwd = os.getcwd()
os.chdir(toolbox)
from general_tools import *
from gaussian_kde_weighted import *
os.chdir(cwd)

#read_dir = '/Users/fiercenator/stuff/projects/02_aging_timescale/gcm_data/CAM_monthly/'

ncfiles = []
ncfiles.extend([read_dir + x for x in os.listdir(read_dir) if x.endswith('.nc')])

monolayer_thicknesses = [1.,2.,4.,8]
m = Basemap(projection='kav7',lon_0=0)

monolayer_thickness = 1.
month = 1.
idx, = np.where([(ncfile.find('aging_L' + str(int(monolayer_thickness)))>0)&(ncfile.find('2010_' + str(int(month)).zfill(2))>0) for ncfile in ncfiles])
ncfile = ncfiles[idx[0]]
#monolayer_thickness = 1.
#idx, = np.where([ncfile.find('aging_L' + str(int(monolayer_thickness)))>0 for ncfile in ncfiles])
#ncfiles12 = [ncfiles[ii] for ii in idx]

f = netCDF4.Dataset(ncfile) 
lon = f.variables['lon'][:]
lat = f.variables['lat'][:]
lev = f.variables['lev'][:]

LON,LAT = np.meshgrid(lon,lat)
con_p = f.variables['con_p'][:]
coag_p = f.variables['coag_p'][:]
tau_p = f.variables['total_p'][:]

con_1 = f.variables['con_m'][:]
coag_1 = f.variables['coag_m'][:]
tau_1 = f.variables['total_m'][:]
mbc_1 = f.variables['bc_pc'][:] + f.variables['bc_accu'][:]


monolayer_thickness = 8.
#idx, =np.where([ncfile.endswith(str(int(monolayer_thickness)) + '.nc') for ncfile in ncfiles])
ncflile = ncfiles[idx[0]]
idx, = np.where([(ncfile.find('aging_L' + str(int(monolayer_thickness)))>0)&(ncfile.find('2010_' + str(int(month)).zfill(2))>0) for ncfile in ncfiles])
ncfile = ncfiles[idx[0]]

f = netCDF4.Dataset(ncfile) 
lon = f.variables['lon'][:]

con_8 = f.variables['con_m'][:]
coag_8 = f.variables['coag_m'][:]
tau_8 = f.variables['total_m'][:]
mbc_8 = f.variables['bc_pc'][:] + f.variables['bc_accu'][:]

fig = plt.figure()
fig.set_size_inches(8.,4.)
ax1 = plt.subplot(2,4,1)
#pcol1 = m.pcolor(LON, LAT, np.log10(np.sum(1/tau_1,axis=0)),latlon=True,clim=[0.,np.log10(24.*7.)],ax=ax1); 
pcol1 = m.pcolor(LON, LAT, np.log10(1/tau_1[55,:,:]),latlon=True,clim=[0.,np.log10(24.*7.)],ax=ax1); 

m.drawcoastlines();
pcol1.set_clim(vmin=np.log10(1/(30.*24.)),vmax=0)

ax2 = plt.subplot(2,4,2)
bbox_vals = ax2.get_position()
ax2.set_position([bbox_vals.x0-0.025,bbox_vals.y0,bbox_vals.width,bbox_vals.height])

#ax3 = plt.subplot(2,3,3)
pcol2 = m.pcolor(LON, LAT, np.log10(1/tau_8[55,:,:]),latlon=True,clim=[0.,np.log10(24.*7.)],ax=ax2); 
#pcol2 = m.pcolor(LON, LAT, np.log10(np.sum(1/tau_8,axis=0)),latlon=True,clim=[0.,np.log10(24.*7.)],ax=ax2); 
m.drawcoastlines();
pcol2.set_clim(vmin=np.log10(1/(30.*24.)),vmax=np.log10(1./6.))
cbar = m.colorbar(pcol2,location='right',ticks=np.log10(1./np.array([24.,24.*7.,24.*30.])))#,ticks=np.log10([1.,24.,24.*7.,24.*30.]))
#cbar.set_ticklabels(['day$^{-1}$','week$^{-1}$','month$^{-1}$'])
cbar.set_ticklabels(['1 day','1 week','1 month'])

ax3 = plt.subplot(2,4,3)
pcol3 = m.pcolor(LON, LAT, np.log10(np.sum(mbc_1,axis=0)),latlon=True,ax=ax3);
pcol3.set_cmap('Greys')
pcol3.set_clim(vmin=-10.5,vmax=-8.5)
#pcol3 = m.pcolor(LON, LAT, np.log10(mbc_8[55,:,:]/mbc_1[55,:,:]),latlon=True,ax=ax3); 
m.drawcoastlines();
cbar3 = m.colorbar(pcol3,location='right')
cbar3.set_ticks([-10,-9,-8])
#cbar3.set_ticklabels([rc('$10^{-10}$'),rc('$10^{-9}$'),rc('$10^{-8}$')])
#cbar3.set_ticklabels(np.array([1e-10,1e-9,1e-8]))
#cbar3.set_ticklabels(10**np.array([-10,-9,-8]))
#cbar4.set_ticks([1.,2.,3.])

ax4 = plt.subplot(2,4,4)
pcol4 = m.pcolor(LON, LAT, np.sum(mbc_8,axis=0)/np.sum(mbc_1,axis=0),latlon=True,clim=[1.,3.],ax=ax4);
pcol4.set_cmap('Purples')
#pcol3 = m.pcolor(LON, LAT, np.log10(mbc_8[55,:,:]/mbc_1[55,:,:]),latlon=True,ax=ax3); 
m.drawcoastlines();
cbar4 = m.colorbar(pcol4,location='right')
pcol4.set_clim(vmin=1.,vmax=4)
#cbar4.set_ticks([1.,2.,3.])

#ax1.set_title(r'$\tau_{\text{age}}$ ' + 'from MAM4,\n 1 monolayer')
#ax2.set_title(r'$\tau_{\text{age}}$ ' + 'from MAM4,\n 8 monolayers')
#ax3.set_title('\n'+ r'$m_{\mathrm{BC},8}/m_{\mathrm{BC},1}$',fontweight='bold')
ax1.set_title(r'$\tau_{\mathrm{age}}$' + '\n1 monolayer')
ax2.set_title(r'$\tau_{\mathrm{age}}$'+ '\n8 monolayers')
ax3.set_title(r'$m_{\mathrm{BC},1}$',fontweight='bold')
ax4.set_title(r'$m_{\mathrm{BC},8}/m_{\mathrm{BC},1}$',fontweight='bold')

bbox_vals = ax3.get_position()
ax3.set_position([bbox_vals.x0+0.03,bbox_vals.y0,bbox_vals.width,bbox_vals.height])

bbox_vals = ax4.get_position()
ax4.set_position([bbox_vals.x0+0.05,bbox_vals.y0,bbox_vals.width,bbox_vals.height])



monolayer_thicknesses = [1.,2.,4.,8]
m = Basemap(projection='kav7',lon_0=0)

monolayer_thickness = 1.
month = 7.
ncflile = ncfiles[idx[0]]
idx, = np.where([(ncfile.find('aging_L' + str(int(monolayer_thickness)))>0)&(ncfile.find('2010_' + str(int(month)).zfill(2))>0) for ncfile in ncfiles])
ncfile = ncfiles[idx[0]]
#monolayer_thickness = 1.
#idx, = np.where([ncfile.find('aging_L' + str(int(monolayer_thickness)))>0 for ncfile in ncfiles])
#ncfiles12 = [ncfiles[ii] for ii in idx]

f = netCDF4.Dataset(ncfile) 
lon = f.variables['lon'][:]
lat = f.variables['lat'][:]
lev = f.variables['lev'][:]

LON,LAT = np.meshgrid(lon,lat)
con_p = f.variables['con_p'][:]
coag_p = f.variables['coag_p'][:]
tau_p = f.variables['total_p'][:]

con_1 = f.variables['con_m'][:]
coag_1 = f.variables['coag_m'][:]
tau_1 = f.variables['total_m'][:]
mbc_1 = f.variables['bc_pc'][:] + f.variables['bc_accu'][:]


monolayer_thickness = 8.
#idx, =np.where([ncfile.endswith(str(int(monolayer_thickness)) + '.nc') for ncfile in ncfiles])
ncflile = ncfiles[idx[0]]
idx, = np.where([(ncfile.find('aging_L' + str(int(monolayer_thickness)))>0)&(ncfile.find('2010_' + str(int(month)).zfill(2))>0) for ncfile in ncfiles])
ncfile = ncfiles[idx[0]]

f = netCDF4.Dataset(ncfile) 
lon = f.variables['lon'][:]

con_8 = f.variables['con_m'][:]
coag_8 = f.variables['coag_m'][:]
tau_8 = f.variables['total_m'][:]
mbc_8 = f.variables['bc_pc'][:] + f.variables['bc_accu'][:]

ax5 = plt.subplot(2,4,5)
#pcol1 = m.pcolor(LON, LAT, np.log10(np.sum(1/tau_1,axis=0)),latlon=True,clim=[0.,np.log10(24.*7.)],ax=ax1); 
pcol5 = m.pcolor(LON, LAT, np.log10(1/tau_1[55,:,:]),latlon=True,clim=[0.,np.log10(24.*7.)],ax=ax5); 

m.drawcoastlines();
pcol5.set_clim(vmin=np.log10(1/(30.*24.)),vmax=0)

ax6 = plt.subplot(2,4,6)
bbox_vals = ax6.get_position()
ax6.set_position([bbox_vals.x0-0.025,bbox_vals.y0,bbox_vals.width,bbox_vals.height])

#ax3 = plt.subplot(2,3,3)
pcol6 = m.pcolor(LON, LAT, np.log10(1/tau_8[55,:,:]),latlon=True,clim=[0.,np.log10(24.*7.)],ax=ax6); 
#pcol2 = m.pcolor(LON, LAT, np.log10(np.sum(1/tau_8,axis=0)),latlon=True,clim=[0.,np.log10(24.*7.)],ax=ax2); 
m.drawcoastlines();
pcol6.set_clim(vmin=np.log10(1/(30.*24.)),vmax=np.log10(1./6.))
cbar = m.colorbar(pcol2,location='right',ticks=np.log10(1./np.array([24.,24.*7.,24.*30.])))#,ticks=np.log10([1.,24.,24.*7.,24.*30.]))
#cbar.set_ticklabels(['day$^{-1}$','week$^{-1}$','month$^{-1}$'])
cbar.set_ticklabels(['1 day','1 week','1 month'])

ax7 = plt.subplot(2,4,7)
pcol7 = m.pcolor(LON, LAT, np.log10(np.sum(mbc_1,axis=0)),latlon=True,ax=ax7);
pcol7.set_cmap('Greys')
#pcol7.set_clim(vmin=-10,vmax=-8)
pcol7.set_clim(vmin=-10.5,vmax=-8.5)
#pcol3 = m.pcolor(LON, LAT, np.log10(mbc_8[55,:,:]/mbc_1[55,:,:]),latlon=True,ax=ax3); 
m.drawcoastlines();
cbar7 = m.colorbar(pcol7,location='right')
cbar7.set_ticks([-10,-9,-8])
#cbar7.set_ticklabels([rc('$10^{-10}$'),rc('$10^{-9}$'),rc('$10^{-8}$')])
#cbar3.set_ticklabels(np.array([1e-10,1e-9,1e-8]))
#cbar3.set_ticklabels(10**np.array([-10,-9,-8]))
#cbar4.set_ticks([1.,2.,3.])

ax8 = plt.subplot(2,4,8)
pcol8 = m.pcolor(LON, LAT, np.sum(mbc_8,axis=0)/np.sum(mbc_1,axis=0),latlon=True,clim=[1.,3.],ax=ax8);
pcol8.set_cmap('Purples')
#pcol3 = m.pcolor(LON, LAT, np.log10(mbc_8[55,:,:]/mbc_1[55,:,:]),latlon=True,ax=ax3); 
m.drawcoastlines();
cbar8 = m.colorbar(pcol8,location='right')
pcol8.set_clim(vmin=1.,vmax=4)
#cbar4.set_ticks([1.,2.,3.])

#ax1.set_title(r'$\tau_{\text{age}}$ ' + 'from MAM4,\n 1 monolayer')
#ax2.set_title(r'$\tau_{\text{age}}$ ' + 'from MAM4,\n 8 monolayers')
#ax3.set_title('\n'+ r'$m_{\mathrm{BC},8}/m_{\mathrm{BC},1}$',fontweight='bold')
#ax5.set_title(r'$\tau_{\mathrm{age}}$' + '\n1 monolayer')
#ax6.set_title(r'$\tau_{\mathrm{age}}$'+ '\n8 monolayers')
#ax7.set_title(r'$m_{\mathrm{BC},1}$',fontweight='bold')
#ax8.set_title(r'$m_{\mathrm{BC},8}/m_{\mathrm{BC},1}$',fontweight='bold')

bbox_vals = ax7.get_position()
ax7.set_position([bbox_vals.x0+0.03,bbox_vals.y0,bbox_vals.width,bbox_vals.height])

bbox_vals = ax8.get_position()
ax8.set_position([bbox_vals.x0+0.05,bbox_vals.y0,bbox_vals.width,bbox_vals.height])

#ax3.set_position([0.55,0.430400989587,0.730913631969,0.574599010413])
#stats_table = table(ax3,
#ax3.set_position([0.55,0.2,0.15,0.6])
fig.savefig('figs/fig3_sens-to-monolayer_month.png',dpi=1000,bbox_inches='tight')