#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 16:30:44 2019

@author: fiercenator
"""

import numpy as np
import os
import netCDF4
import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap

# toolbox = '../../code/process_partmc/'
toolbox = '../../../process_partmc/'

cwd = os.getcwd()
os.chdir(toolbox)
from general_tools import *
from gaussian_kde_weighted import *
os.chdir(cwd)

#read_dir = '/Users/fiercenator/stuff/projects/02_aging_timescale/gcm_data/CAM_aging/'
read_dir = '/Users/fier887/Downloads/02_aging_timescale/gcm_data/CAM_aging/'

ncfiles = []
cwd = os.getcwd()
os.chdir(read_dir)
ncfiles.extend([read_dir + x for x in os.listdir(read_dir) if x.endswith('.nc')])
os.chdir(cwd)

#monolayer_thicknesses = [1.,2.,4.,8]
#m = Basemap(projection='kav7',lon_0=0)

monolayer_thicknesses = [8.,4.,2.,1.]#[1.,2.,4.,8.]#
mm = 1

monolayer_thickness = 1.
idx, =np.where([ncfile.endswith(str(int(monolayer_thickness)) + '.nc') for ncfile in ncfiles])
ncfile = ncfiles[idx[0]]

f = netCDF4.Dataset(ncfile) 
lon = f.variables['lon'][:]
lat = f.variables['lat'][:]


LON,LAT = np.meshgrid(lon,lat)


fig = plt.figure()
gs = plt.GridSpec(2,1,height_ratios = [2,1],hspace=0.5)
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])
#ax3 = plt.subplot(gs[2])
#hln1 = list(np.zeros(4))
for monolayer_thickness in monolayer_thicknesses:
    idx, =np.where([ncfile.endswith(str(int(monolayer_thickness)) + '.nc') for ncfile in ncfiles])
    ncfile = ncfiles[idx[0]]
    
    f = netCDF4.Dataset(ncfile)
    tau_m = f.variables['total_m'][:]
    tau_p = f.variables['total_p'][:]
    mbc_fresh = f.variables['bc_pc'][:]
    
    LEV = np.zeros(mbc_fresh.shape);
    ii = 0
    for lev in f.variables['lev']:
        LEV[ii,:,:] = lev
        ii = ii + 1
    idx1,idx2,idx3 = np.where((mbc_fresh[:]>0.)&(LEV>200))
    #idx1,idx2,idx3 = np.where(mbc_fresh[:]>0.)    
    aging_rate_diff =(1./tau_m[idx1,idx2,idx3] - 1./tau_p[idx1,idx2,idx3])#/(1./tau_p[:]) #np.reshape((1./tau_m[:] - 1./tau_p[:])/(1./tau_p[:]),np.prod(np.shape(tau_p[:]))).data

#    print(sum(mbc_fresh[idx1,idx2,idx3]*abs(aging_rate_diff)/(1./tau_p[idx1,idx2,idx3]))/sum(mbc_fresh[idx1,idx2,idx3]))
#    print(sum(mbc_fresh[idx1,idx2,idx3]*(aging_rate_diff)/(1./tau_p[idx1,idx2,idx3]))/sum(mbc_fresh[idx1,idx2,idx3]))    
#    print(sum(mbc_fresh[idx1,idx2,idx3].data*abs(aging_rate_diff)**2)/sum(mbc_fresh[idx1,idx2,idx3].data))
    
    kernel_1d = gaussian_kde_weighted(aging_rate_diff/(1./tau_p[idx1,idx2,idx3]),weights=mbc_fresh[idx1,idx2,idx3],bw_method=0.1)
    hln1, = ax1.plot(np.linspace(-1.,6.,200),kernel_1d(np.linspace(-1.,6.,200)))
#    kernel_1d = gaussian_kde_weighted(abs(aging_rate_diff)/(1./tau_p[idx1,idx2,idx3]),weights=mbc_fresh[idx1,idx2,idx3],bw_method=0.1)
#    hln2, = ax2.plot(np.linspace(-1.,6.,200),kernel_1d(np.linspace(-1.,6.,200)),color=hln1.get_color())

#    ax2.barh(np.log10(monolayer_thickness),sum(mbc_fresh[idx1,idx2,idx3]*abs(aging_rate_diff)/(1./tau_p[idx1,idx2,idx3]))/sum(mbc_fresh[idx1,idx2,idx3]),color='white',edgecolor=hln1.get_color(),height=0.2)
    ax2.barh(np.log10(monolayer_thickness),sum(mbc_fresh[idx1,idx2,idx3]*(aging_rate_diff)/(1./tau_p[idx1,idx2,idx3]))/sum(mbc_fresh[idx1,idx2,idx3]),color=hln1.get_color(),height=0.2)

#    ax3.scatter(sum(mbc_fresh[idx1,idx2,idx3]*abs(aging_rate_diff)/(1./tau_p[idx1,idx2,idx3]))/sum(mbc_fresh[idx1,idx2,idx3]),monolayer_thickness,edgecolor=hln1.get_color(),marker='o',facecolor='w')
#    ax2.scatter(sum(mbc_fresh[idx1,idx2,idx3]*(aging_rate_diff)/(1./tau_p[idx1,idx2,idx3]))/sum(mbc_fresh[idx1,idx2,idx3]),monolayer_thickness,edgecolor=hln1.get_color(),marker='o',facecolor=hln1.get_color())  
#    ax1.set_xlim([-1,4])
#    ax2.set_xlim([-1,4])    
#    ax3.set_xlim([-1,4])    
#ax2.set_yscale('log')

ylims = ax1.get_ylim(); ylims= (0.,ylims[1])
ax1.plot([0.,0.],ylims,color='k',linewidth=0.5);
ax1.set_ylim(ylims)
ylims = ax2.get_ylim()
ax2.plot([0.,0.],ylims,color='k',linewidth=0.5);
ax2.set_ylim(ylims)
#ax3.plot([0.,0.],ylims,color='k',linewidth=0.5);
#ax3.set_ylim(ylims)

#ax1.xaxis.tick_top()
#ax2.xaxis.tick_top()
#ax2.set_xticklabels('')
#ax2.spines['right'].set_visible(False)
#ax2.spines['bottom'].set_visible(False)
#ax3.spines['right'].set_visible(False)
#ax3.spines['top'].set_visible(False)

ax1.set_xlim([-1,4])
ax2.set_xlim(ax1.get_xlim())
#ax3.set_xlim(ax1.get_xlim())

ax1.set_xticklabels([str(tick*100) + '%' for tick in ax1.get_xticks()])
ax2.set_xticklabels([str(tick*100) + '%' for tick in ax2.get_xticks()])

ax1.set_xlabel('relative difference between mixing rate and aging rate')
ax2.set_xlabel('mean relative difference between mixing rate and aging rate')
leg_strs = list(range(len(monolayer_thicknesses)))
for mm in range(len(monolayer_thicknesses)):
    if monolayer_thicknesses[mm] == 1:
        leg_strs[mm] = str(int(monolayer_thicknesses[mm])) + ' monolayer'
    else:
        leg_strs[mm] = str(int(monolayer_thicknesses[mm])) + ' monolayers'
ax1.legend(leg_strs)
ax2.get_yaxis().set_visible(False)
#ax2.spines['top'].set_visible(False)

#ax1.legend(['1 monolayer', '2 monolayers', '4 monolayers', '8 monolayers'])

fig.savefig('../../figs/fig5_optimize-monolayer.png',dpi=1000)

# fig.savefig('figs/fig5_optimize-monolayer.pdf',dpi=1000)
#fig.savefig('figs/fig5_optimize-monolayer_with-mean-abs-err.png')

