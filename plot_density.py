#!/usr/bin/env python2
import sys
import linecache

import numpy as np
import os

import matplotlib as mpl
mpl.use('Agg')   # to avoid window

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

AU = 1.496e13

filename = "dust_density.binp"
print filename

#################################
############### read the grids

grid = linecache.getline('amr_grid.inp', 6)
print grid
ngrids = [int(s) for s in grid.split() if s.isdigit()]
nrad  = ngrids[0]
nthet = ngrids[1]
nphi  = ngrids[2]
print "nrad", nrad
print "nthet", nthet
print "nphi", nphi

grids = np.genfromtxt('amr_grid.inp', skip_header=6, dtype=float)
print grids.shape

nrad_b  = nrad+1
nthet_b = nthet+1
nphi_b  = nphi+1
print "ngrids[0]", ngrids[0]
print "ngrids[1]", ngrids[1]
print "ngrids[2]", ngrids[2]

nrad_b1 = grids[0:nrad:1]
nrad_b2 = grids[1:nrad+1:1]
rad = (nrad_b1+nrad_b2)/2.0/AU

nthet_b1 = grids[nrad+1:nrad+1+nthet:1]
nthet_b2 = grids[nrad+2:nrad+2+nthet:1]
thet = (nthet_b1+nthet_b2)/2.0

nphi_b1 = grids[nrad+2+nthet:nrad+2+nthet+nphi:1]
nphi_b2 = grids[nrad+3+nthet:nrad+3+nthet+nphi:1]
phi = (nphi_b1+nphi_b2)/2.0


#############################################
############### read the density data

head = np.fromfile(filename, count=4, dtype=int)
print "head[0]", head[0]
print "head[1]", head[1]
print "head[2]", head[2]
print "head[3]", head[3]

if head[2]!=(nrad*nthet*nphi):
    print ' ERROR'
    print ' Number of grid points in '+filename+' is different from that in amr_grid.inp'
    quit()

print " Reading the dust_density file"

if head[1] == 8:
    f = open(filename, "rb")  # reopen the file
    f.seek(32, os.SEEK_SET)  # seek
    density_3d = np.fromfile(f, count=-1, dtype=np.float64)
elif head[1]==4:
    f = open(filename, "rb")  # reopen the file
    f.seek(32, os.SEEK_SET)  # seek
    density_3d = np.fromfile(f, count=-1, dtype=np.float32)
else:
    print 'ERROR'
    print 'Unknown datatype in '+filename
f.close()

print " dust_density reading: done."

density_3d = np.array(density_3d).reshape(nphi,nthet,nrad)
print " density_3d transformed to nphi,nthet,nrad"


################3333333333333333333333
##### read the grid

"""

surf_dens2d = np.zeros((nthet,nrad), dtype=float)

const=np.sin(dphi)*1.496e13

for i in range(nthet):
    print i
    for k in range(nrad):
        for j in range(nphi):
            surf_dens2d[i,k]=surf_dens2d[i,k] + density_3d[j,i,k]*rad[k]*const

"""

plt.pcolormesh(rad, phi, density_3d[:,0,:]) #, norm=LogNorm())
plt.xlabel('Radial Distance [AU]')
plt.ylabel('Polar Angle')
#plt.colorbar()
#r.set_label('Vertical Density of Dust at phi = 0 [ g/cm^2 ]', rotation=90)
plt.savefig('verticaldens.png')

quit()



