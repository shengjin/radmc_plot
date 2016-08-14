#!/usr/bin/env python2
import sys
import linecache

import numpy as np
import os

import matplotlib as mpl
mpl.use('Agg')   # to avoid the window

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

#sys.path.append('/turquoise/users/sjin/python/')

debug=True

PIhalf = 1.57079637051

AU = 1.496e13

n = ""
file_temp = "%s%s%s" % ("dust_temperature", n, ".bdat")
file_vgas = "%s%s%s" % ("gas_velocity", n, ".inp")
#file_vgas = "%s%s%s" % ("dust_density_small", n, ".binp")


#################################
############### read the grids

grid = linecache.getline('amr_grid.inp', 6)
print grid
ngrids = [int(s) for s in grid.split() if s.isdigit()]
nrad  = ngrids[0]
nphi = ngrids[1]
nz  = ngrids[2]
print "nrad", nrad
print "nthet", nthet
print "nphi", nphi

grids = np.genfromtxt('amr_grid.inp', skip_header=6, dtype=float)
print grids.shape

nrad_b  = nrad+1
nphi_b = nphi+1
nz_b  = nz+1
print "ngrids[0]", ngrids[0]
print "ngrids[1]", ngrids[1]
print "ngrids[2]", ngrids[2]

nrad_b1 = grids[0:nrad:1]
nrad_b2 = grids[1:nrad+1:1]
rad = (nrad_b1+nrad_b2)/2.0/AU

nphi_b1 = grids[nrad+1:nrad+1+nphi:1]
nphi_b2 = grids[nrad+2:nrad+2+nphi:1]
phi = (nphi_b1+nphi_b2)/2.0

nz_b1 = grids[nrad+2+nphi:nrad+2+nphi+nz:1]
nz_b2 = grids[nrad+3+nphi:nrad+3+nphi+nz:1]
zz = (nz_b1+nz_b2)/2.0


##########################
##### read the gas_velocity.inp

print " Reading the gas_velocity file"
gas_velocity=np.genfromtxt(file_vgas, skip_header=2, dtype=float)
print " gas_velocity reading: done."

m_gas,n_gas = gas_velocity.shape
if m_gas!=(nrad*nphi*nz):
    print ' ERROR'
    print ' Number of grid points in '+file_vgas+' is different from that in amr_grid.inp'
if n_gas!=3:
    print ' ERROR'
    print ' Number of columns in '+file_vgas+' is not equal to 3'
    quit()


vgas_r = np.zeros(m_gas)
vgas_z = np.zeros(m_gas)
vgas_phi = np.zeros(m_gas)
vgas_r = gas_velocity[:,0]
vgas_z = gas_velocity[:,1]
vgas_phi = gas_velocity[:,2]

vgas_r = np.array(vgas_r).reshape(nrad,nphi,nz)
vgas_z = np.array(vgas_z).reshape(nrad,nphi,nz)
vgas_phi = np.array(vgas_phi).reshape(nrad,nphi,nz)


### output the vgas_r and vgas_phi for check.
if debug:
    plt.pcolormesh(rad, phi, vgas_r[:,:,0].T, vmin=vgas_r[:,:,0].min(), vmax=vgas_r[:,:,0].max())
    plt.colorbar()
    plt.savefig('vgas-r_r_phi.png')
    plt.clf()
    plt.pcolormesh(rad, z, vgas_r[:,0,:].T, vmin=vgas_r[:,0,:].min(), vmax=vgas_r[:,0,:].max())
    plt.colorbar()
    plt.savefig('vgas-r_r_z.png')
    plt.clf()
    plt.pcolormesh(rad, z, vgas_phi[:,0,:].T, vmin=vgas_phi[:,0,:].min(), vmax=vgas_phi[:,0,:].max())
    plt.colorbar()
    plt.savefig('vgas-phi_r_z.png')
    plt.clf()
    plt.pcolormesh(rad, phi, vgas_phi[:,:,0].T, vmin=vgas_phi[:,:,0].min(), vmax=vgas_phi[:,:,0].max())
    plt.colorbar()
    plt.savefig('vgas-phi_r_phi.png')
    plt.clf()

