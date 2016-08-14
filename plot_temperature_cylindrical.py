#!/usr/bin/env python2
import sys
import linecache

import numpy as np
import os

import matplotlib as mpl
mpl.use('Agg')   # to avoid the window

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

sys.path.append('/turquoise/users/sjin/python/')

PIhalf = 1.57079637051

AU = 1.496e13

n = ""
filename = "%s%s%s" % ("dust_temperature", n, ".bdat")
print filename

#################################
############### read the grids

grid = linecache.getline('amr_grid.inp', 6)
print grid
ngrids = [int(s) for s in grid.split() if s.isdigit()]
nrad  = ngrids[0]
nphi = ngrids[1]
nz  = ngrids[2]
print "nrad", nrad
print "nphi", nphi
print "nz", nz

grids = np.genfromtxt('amr_grid.inp', skip_header=6, dtype=float)
print grids.shape

nrad_b1 = grids[0:nrad:1]
nrad_b2 = grids[1:nrad+1:1]
rad = (nrad_b1+nrad_b2)/2.0/AU

nphi_b1 = grids[nrad+1:nrad+1+nphi:1]
nphi_b2 = grids[nrad+2:nrad+2+nphi:1]
phi = (nphi_b1+nphi_b2)/2.0

nz_b1 = grids[nrad+2+nphi:nrad+2+nphi+nz:1]
nz_b2 = grids[nrad+3+nphi:nrad+3+nphi+nz:1]
zz = (nz_b1+nz_b2)/2.0


#############################################
############### read the density data

head = np.fromfile(filename, count=4, dtype=int)
print "head[0]", head[0]
print "head[1]", head[1]
print "head[2]", head[2]
print "head[3]", head[3]


print " Reading the dust_temp file"

if head[1] == 8:
    f = open(filename, "rb")  # reopen the file
    f.seek(32, os.SEEK_SET)  # seek
    temp_3d = np.fromfile(f, count=-1, dtype=np.float64)
elif head[1]==4:
    f = open(filename, "rb")  # reopen the file
    f.seek(32, os.SEEK_SET)  # seek
    temp_3d = np.fromfile(f, count=-1, dtype=np.float32)
else:
    print 'ERROR'
    print 'Unknown datatype in '+filename
f.close()

print " dust_temp reading: done."

temp_3d = np.array(temp_3d).reshape(nrad,nphi,nz)


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

plt.pcolormesh(rad, zz, temp_3d[:,0,:].T)
#plt.ylim([1.6,0.9])
#plt.pcolormesh(rad, phi, temp_3d[:,0,:])
#plt.pcolormesh(rad, thet, temp_3d[:,nthet-30,:], vmin=30, vmax=80)
#plt.pcolormesh(rad, phi, temp_3d[:,nthet-1,:], vmin=temp_3d[0,:,:].min(), vmax=temp_3d[0,:,:].max())
#plt.xlabel('Radial Distance [AU]')
#plt.ylabel('Phi')
plt.colorbar()
#r.set_label('Vertical Density of Dust at phi = 0 [ g/cm^2 ]', rotation=90)
plt.savefig('verticaldens.png')

plt.clf()



###### Plot temp at the radial distance for one slice

onephi=0
thet_nth = nthet - 1 # midplane

oneradial = np.zeros(nrad)
oneradial = temp_3d[onephi,thet_nth,:] 
np.savetxt('au_T.out', np.transpose([rad,oneradial]))

#plt.ylim([5,150])
#thet = "%0.1f" % np.degrees(thet[thet_nth])
label = str(thet)
label='T_dust (phi=0, thet=%s )' % label
#plt.loglog(rad,oneradial, label="mid-plane", color="red")

print thet
print temp_3d[onephi,:,5]
plt.plot(thet,temp_3d[onephi,:,5], label='T_dust (surface)')
#plt.plot(rad,temp_3d[onephi,0,:], label='T_dust (surface)')
#plt.loglog(rad,temp_3d[onephi,0,:], label='T_dust (surface)')
#plt.plot(phi,temp_3d[:,0,80], label="mid-plane", color="red")
#plt.plot(phi,temp_3d[:,10,80], label="mid-plane", color="red")
#plt.plot(rad,temp_3d[onephi,thet_nth,:], label="mid-plane", color="red")


##### Plot the temperature in hydro simulation
zeta = 0.43
scaling_radius = 10

Y = 0.24
X = 1.0-Y
Mu = 1.0/(X/2.006+Y/4.008)
Mp = 1.67262158E-24
kb = 1.380658E-16
#temp = cs_square*Mu*MP/kb
H0 = 0.05
r0 = 10.0*1.4959787E+13
Mstar = 0.55*1.9891E+33
G = 6.67259E-8
vk0_square = G*Mstar/r0

Const = H0**2.0*vk0_square*Mu*Mp/kb

Tgas = np.zeros(nrad)
Tgas = (rad[:]/scaling_radius)**-zeta*Const
Tgas_2 = (rad[:]/scaling_radius)**-0.5*Const
#plt.loglog(rad,Tgas,label='T_gas', color='black', lw=1, ls='-.')
#plt.xlim([4,120])
#plt.plot(rad,Tgas,label='T_gas', color='green', lw=1, ls='-.')
#plt.plot(rad,Tgas_2,label='T_gas_radmc3d', color='black', lw=1)
#plt.plot(rad,Tgas_2*0+75.0,label='T_gas_radmc3d', color='black', lw=1)

plt.legend()
plt.xlabel('Semi-major Axis ( AU )')
plt.ylabel('Temperature ( K ) ')

plt.savefig('temp_slice.png')

