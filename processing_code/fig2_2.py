#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 15:28:21 2019

@author: fiercenator

Figure 2: timescale parameterization in GCM
"""


import numpy as np
import os
import netCDF4
import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
from basemap import Basemap
from matplotlib.colors import LogNorm
from scipy.signal import savgol_filter
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker
from gaussian_kde_weighted import *

def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, np.sqrt(variance))

read_dir = '/Users/fier887/Downloads/02_aging_timescale/gcm_data/CAM_aging/'
toolbox = '../../../process_partmc/'
#toolbox = '../../code/process_partmc/'
cwd = os.getcwd()
os.chdir(toolbox)
from general_tools import *
from gaussian_kde_weighted import *
os.chdir(cwd)

#read_dir = '/Users/fiercenator/stuff/projects/02_aging_timescale/gcm_data/CAM_aging/'

ncfiles = []
#cwd = os.getcwd()
#os.chdir(read_dir)
ncfiles.extend([read_dir + x for x in os.listdir(read_dir) if x.endswith('.nc')])
#os.chdir(cwd)

monolayer_thicknesses = [1.,2.,4.,8]
m = Basemap(projection='kav7',lon_0=0)


monolayer_thickness = 1.
idx, =np.where([ncfile.endswith(str(int(monolayer_thickness)) + '.nc') for ncfile in ncfiles])
ncfile = ncfiles[idx[0]]

f = netCDF4.Dataset(ncfile) 
lon = f.variables['lon'][:]
lat = f.variables['lat'][:]


LON,LAT = np.meshgrid(lon,lat)
mbc_fresh = f.variables['bc_pc'][:]
con_p = f.variables['con_p'][:]
coag_p = f.variables['coag_p'][:]
tau_p = f.variables['total_p'][:]

fig = plt.figure()
ax1 = plt.subplot(2,2,1)
ax3 = plt.subplot(2,2,2)
ax2 = plt.subplot(2,2,3)
ax4 = plt.subplot(2,2,4)
#gs1 = plt.GridSpec(5,2)
#gs1.update(hspace=-0.4)
#ax1 = plt.subplot(gs1[:2,:])
#ax2 = plt.subplot(gs1[4,0])
#ax3 = plt.subplot(gs1[4,1])
#ax1 = plt.subplot(1,3,1)f

pcol1 = m.pcolor(LON, LAT, np.log10(1./tau_p[55,:,:]),latlon=True,clim=[0.,np.log10(24.*7.)],ax=ax1); 
m.drawcoastlines(ax=ax1);
pcol1.set_clim(vmin=np.log10(1./(30.*24.)),vmax=np.log10(1./6.))

ax1.set_title(r'$\tau_{\mathrm{mix}}$',fontweight='bold')

cbar1 = m.colorbar(pcol1,location='right',pad="10%",ticks=np.log10(1./np.array([24.,24.*7.,24.*30.])),ax=ax1)
cbar1.set_ticklabels(['1 day','1 week','1 month'])

#fig.savefig('figs/fig2_timescales_overall.png')

#fig = plt.figure()
#ax2 = plt.subplot(1,3,2)
#pcol2 = m.pcolor(LON, LAT, np.log10(1/con_p[55,:,:]),latlon=True,clim=[0.,1.],ax=ax2); 
pcol2 = m.pcolor(LON, LAT, (1./con_p[55,:,:])/(1./tau_p[55,:,:]),latlon=True,clim=[0.,1.],ax=ax2); 
m.drawcoastlines(ax=ax2);
pcol2.set_cmap('Blues')
pcol2.set_clim(vmin=0.,vmax=1.)
#pcol2.set_clim(vmin=np.log10(1./(30.*24.)),vmax=np.log10(1./6.))

##ax3 = plt.subplot(1,3,3)
#pcol3 = m.pcolor(LON, LAT, (1./coag_p[55,:,:])/(1./tau_p[55,:,:]),latlon=True,clim=[0.,1.],ax=ax3); 
##pcol3 = m.pcolor(LON, LAT, np.log10(1/coag_p[55,:,:]),latlon=True,clim=[0.,np.log10(24.*7.)],ax=ax3); 
#m.drawcoastlines(ax=ax3);
##pcol3.set_clim(vmin=np.log10(1./(30.*24.)),vmax=np.log10(1./6.))
#pcol3.set_cmap('Reds')
#pcol3.set_clim(vmin=0.,vmax=1.)

ax2.set_title(r'$f_{\mathrm{cond}}$')
#ax2.set_title('condensation\'s contribution')

#ax3.set_ylabel('frequency dist. of BC mass')
ax3.set_ylabel('frequency dist., weighted\nby mixing ratio of fresh BC')

#ax4.set_ylabel(r'$f_{\mathrm{cond}}$')
#ax3.set_title('coagulation\'s\ncontribution')
#ax2.set_title(r'$1/\tau_{\mathrm{cond}}$')
#ax3.set_title(r'$1/\tau_{\mathrm{coag}}$')

cbar2 = m.colorbar(pcol2,location='right',pad="10%",ax=ax2,ticks=(0.,0.5,1.))
cbar2.set_ticklabels(['0%','50%','100%'])

LEV = np.zeros(mbc_fresh.shape);
ii = 0
for lev in f.variables['lev']:
    LEV[ii,:,:] = lev
    ii = ii + 1

idx1,idx2,idx3 = np.where((mbc_fresh[:]>0.)&(LEV>200))
f_cond = tau_p[idx1,idx2,idx3]/con_p[idx1,idx2,idx3].ravel()
f_coag = tau_p[idx1,idx2,idx3]/coag_p[idx1,idx2,idx3].ravel()
mbc = mbc_fresh[idx1,idx2,idx3]
tau_bins = np.logspace(-1,3,80)
contribution_bins = np.linspace(0.,1.,40)

idx = np.digitize(tau_p[idx1,idx2,idx3],tau_bins)
hist1d = np.zeros(len(tau_bins))
fcond_mean = np.zeros(len(tau_bins))
fcond_std = np.zeros(len(tau_bins))
fcoag_mean = np.zeros(len(tau_bins))
fcoag_std = np.zeros(len(tau_bins))
dtau = np.log10(tau_bins[1]) - np.log10(tau_bins[0])
for ii in range(len(tau_bins)):
    hist1d[ii] = sum(mbc[idx==ii])/(sum(mbc)*dtau)
    if hist1d[ii]>0:
        fcond_mean[ii],fcond_std[ii] = weighted_avg_and_std(f_cond[idx==ii], mbc[idx==ii])
        fcoag_mean[ii],fcoag_std[ii] = weighted_avg_and_std(f_coag[idx==ii], mbc[idx==ii])

kernel_1d = gaussian_kde_weighted(np.log10(tau_p[idx1,idx2,idx3]),weights=mbc_fresh[idx1,idx2,idx3],bw_method=0.1)
ax3.plot(np.log10(tau_bins),kernel_1d(np.log10(tau_bins)),color='k')

col_cond = '#0033cc'
col_coag = '#009999'

ax4.errorbar(np.log10(tau_bins),savgol_filter(fcond_mean,5,1),savgol_filter(fcond_std,5,1),color=col_cond); 
ax4.errorbar(np.log10(tau_bins),savgol_filter(fcoag_mean,5,1),savgol_filter(fcoag_std,5,1),color=col_coag); 
ax4.set_xticks(np.log10([1,24,24*7,24*30])); 
ax4.set_xticklabels(['1 hour', '1 day', '1 week', '1 month'],rotation=45)
ax4.set_xlabel(r'$\tau_{\mathrm{mix}}$')

ax3.set_xticks(np.log10([1,24,24*7,24*30])); 
ax3.set_xticklabels([])

ybox1 = TextArea(r'$f_{\mathrm{cond}}$', textprops=dict(color=col_cond, rotation=90,ha='left',va='bottom'))
ybox2 = TextArea(', ', textprops=dict(color='k', rotation=90,ha='left',va='bottom'))
ybox3 = TextArea(r'$f_{\mathrm{coag}}$', textprops=dict(color=col_coag, rotation=90,ha='left',va='bottom'))
ybox = VPacker(children=[ybox1, ybox2, ybox3],align="bottom", pad=0, sep=0)
anchored_ybox = AnchoredOffsetbox(loc=8, child=ybox, pad=0., frameon=False, bbox_to_anchor=(-0.295, 0.09), 
                                  bbox_transform=ax4.transAxes, borderpad=0.)
ax4.add_artist(anchored_ybox)

 
ax4.set_xlim(np.log10([0.5,24*30*1.01]))
ax3.set_xlim(np.log10([0.5,24*30*1.01]))

ax3.set_position([0.75,0.47,0.3,0.36])
ax4.set_position([0.75,0.16,0.3,0.26])

fig.savefig('figs/fig2_combined.pdf',dpi=1000,bbox_inches='tight')

# fig.savefig('figs/fig2_combined.png',dpi=1000,bbox_inches='tight')
#ax3.set_position([0.7,0.55,0.3,0.3])
#ax4.set_position([0.7,0.15,0.3,0.28])


#kde_cond = gaussian_kde_weighted([np.log10(tau_p[idx1,idx2,idx3]).ravel(),tau_p[idx1,idx2,idx3].ravel()/con_p[idx1,idx2,idx3].ravel()],bw_method=0.03,weights=mbc_fresh[idx1,idx2,idx3].ravel())
#kde_coag = gaussian_kde_weighted([np.log10(tau_p[idx1,idx2,idx3]).ravel(),tau_p[idx1,idx2,idx3].ravel()/coag_p[idx1,idx2,idx3].ravel()],weights=mbc_fresh[idx1,idx2,idx3].ravel())
#kde_1d = gaussian_kde_weighted(np.log10(tau_p[idx1,idx2,idx3]).ravel(),weights=mbc_fresh[idx1,idx2,idx3].ravel())
#hist1d = np.zeros(len(tau_bins))
#hist2d_cond = np.zeros([len(contribution_bins),len(tau_bins)])
#hist2d_coag = np.zeros([len(contribution_bins),len(tau_bins)])
#for ii in range(len(tau_bins)):
#    hist1d[ii] = kde_1d(tau_bins[ii])
#    for jj in range(len(contribution_bins)):
#        hist2d_cond[jj,ii] = kde_cond([contribution_bins[jj],np.log10(tau_bins[ii])])
#        hist2d_coag[jj,ii] = kde_coag([contribution_bins[jj],np.log10(tau_bins[ii])])
    
#kernel_1d = gaussian_kde_weighted(np.log10(1./tau_p[idx1,idx2,idx3]),weights=mbc_fresh[idx1,idx2,idx3],bw_method=0.1)
#kernel_2d = gaussian_kde_weighted(np.vstack([np.log10(1./tau_p[idx1,idx2,idx3]),tau_p[idx1,idx2,idx3]/con_p[idx1,idx2,idx3]]),weights=mbc_fresh[idx1,idx2,idx3],bw_method=0.1)

#dat = ax3.hist(np.log10(tau_p[idx1,idx2,idx3]),bins=100,weights=mbc_fresh[idx1,idx2,idx3])
#hist1d = dat[0]
#log10tau_edges= dat[1]
#dat = ax4.hist2d(np.log10(tau_p[idx1,idx2,idx3]),tau_p[idx1,idx2,idx3]/con_p[idx1,idx2,idx3],bins=log10tau_edges,weights=mbc_fresh[idx1,idx2,idx3])
#hist2d_cond = dat[0]
#contribution_edges = dat[2]
#contribution_mids = (contribution_edges[:-1] + contribution_edges[1:])/2
#fcond_lb = np.zeros(np.shape(hist1d))
#fcond_ub = np.zeros(np.shape(hist1d))
#
#for ii in range(hist2d_cond.shape[0]):
#    hist2d_cond[ii,:] = hist2d_cond[ii,:]/hist1d[ii]
#    
#    thresh = 0.25;
#    idx, = np.where(np.cumsum(hist2d_cond[ii,:])<=thresh)
#    idx_these = idx[-1]
#    idx, = np.where(np.cumsum(hist2d_cond[ii,:])>=thresh)
#    idx_these = [idx_these,idx[0]]
#    fcond_lb[ii] = np.interp(thresh, np.cumsum(hist2d_cond[ii,:]), contribution_mids)
#    
#    thresh = 0.75;
#    idx, = np.where(np.cumsum(hist2d_cond[ii,:])<=thresh)
#    idx_these = idx[-1]
#    idx, = np.where(np.cumsum(hist2d_cond[ii,:])>=thresh)
#    idx_these = [idx_these,idx[0]]
#    fcond_ub[ii] = np.interp(thresh, np.cumsum(hist2d_cond[ii,:]), contribution_mids)
#    
#dat = ax4.hist2d(np.log10(tau_p[idx1,idx2,idx3]),tau_p[idx1,idx2,idx3]/coag_p[idx1,idx2,idx3],bins=log10tau_edges,weights=mbc_fresh[idx1,idx2,idx3])
##dat = ax4.hist2d(np.log10(tau_p[idx1,idx2,idx3]),np.log10(coag_p[idx1,idx2,idx3]),bins=log10tau_edges,weights=mbc_fresh[idx1,idx2,idx3])
#hist2d_coag = dat[0]
#
#
#fcoag_lb = np.zeros(np.shape(hist1d))
#fcoag_ub = np.zeros(np.shape(hist1d))
#for ii in range(hist2d_coag.shape[0]):
#    hist2d_coag[ii,:] = hist2d_coag[ii,:]/hist1d[ii]
#    
#    thresh = 0.25;
#    idx_these = [idx_these,idx[0]]
#    fcoag_lb[ii] = np.interp(thresh, np.cumsum(hist2d_coag[ii,:]), contribution_mids)
#    
#    thresh = 0.75;
#    fcoag_ub[ii] = np.interp(thresh, np.cumsum(hist2d_coag[ii,:]), contribution_mids)

    
#ax3.hist(np.log10(1./tau_p.ravel()),bins=100,weights=mbc_fresh.ravel())
#h2d = ax4.hist2d(np.log10(1./tau_p[idx1,idx2,idx3]),tau_p[idx1,idx2,idx3]/con_p[idx1,idx2,idx3],weights=mbc_fresh[idx1,idx2,idx3],bins=100,vmin=0,vmax=1e-9);
#ax4.pcolor(h2d[1],h2d[2],h2d[0],norm=LogNorm(vmin=h2d[0][h2d[0]>0].min(), vmax=h2d[0].max()))
#ax4.scatter(np.log10(tau_p[idx1,idx2,idx3]),np.log10(tau_p[idx1,idx2,idx3]),0.1,'k')

#ax4.plot(np.linspace(-2,4,100),np.linspace(-2,4,100),color='k',linestyle='--')
#ax4.scatter(np.log10(tau_p[idx1,idx2,idx3]),np.log10(coag_p[idx1,idx2,idx3]),0.01,alpha=0.01)
#ax4.scatter(np.log10(tau_p[idx1,idx2,idx3]),np.log10(con_p[idx1,idx2,idx3]),0.01,alpha=0.01)
#ax4.hist(np.log10(coag_p[idx1,idx2,idx3]),bins=100,weights=mbc_fresh[idx1,idx2,idx3])
# ax4.hist(con_p[idx1,idx2,idx3],bins=100,weights=mbc_fresh[idx1,idx2,idx3])
#ax4.scatter(tau_p[idx1,idx2,idx3],tau_p[idx1,idx2,idx3]/con_p[idx1,idx2,idx3],1); plt.xscale('log');
#ax4.scatter(tau_p[idx1,idx2,idx3],tau_p[idx1,idx2,idx3]/con_p[idx1,idx2,idx3],1); plt.xscale('log');

#ax3.hist(np.log10(1./tau_p[idx1,idx2,idx3]),bins=100,weights=mbc_fresh[idx1,idx2,idx3],weights=mbc_fresh[idx1,idx2,idx3])
#ax4.hist(np.log10(1./tau_p[idx1,idx2,idx3]),tau_p[idx1,idx2,idx3]/con_p[idx1,idx2,idx3],weights=mbc_fresh[idx1,idx2,idx3]);
#fig.savefig('figs/fig2_timescales.png',dpi=1000)
#fig.savefig('figs/fig2_timescales_contribution.png')
#fig.savefig('figs/fig2_timescales_cb1.png')


#fig.savefig('figs/fig2_timescales_cb2.png')
