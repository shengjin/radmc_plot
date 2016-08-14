#!/usr/bin/env python2
import os
import sys
import linecache
import time
import numpy as np
import math

from matplotlib.colors import LogNorm
import matplotlib as mpl
mpl.use('Agg')   # to avoid the window

#then
import matplotlib.pyplot as plt

#wavelength = 10
#smallorbig = "small"
wavelength = 1000
smallorbig = "large"


grid = linecache.getline('amr_grid.inp', 6)
print grid
ngrids = [int(s) for s in grid.split() if s.isdigit()]
print "ngrids[0]", ngrids[0]
print "ngrids[1]", ngrids[1]
print "ngrids[2]", ngrids[2]
nrad=ngrids[0]
nthet=ngrids[1]
nphi=ngrids[2]
print "nrad", nrad
print "nthet", nthet
print "nphi", nphi


rad_min = 0.24
rad_max = 15.6
scaling_radius = 10.0
rin  = rad_min  * scaling_radius
rout = rad_max  * scaling_radius
#  // radius scale: LIN or LOG
scale  = "LIN" 
#  // define the radial grid
rad = np.zeros(nrad)
if (scale == "LIN"):
    dr = (rout-rin)/(nrad-1) 
    for i in range(nrad):
        rad[i] = rin  + i*dr
if (scale == "LOG"):
    dr = (rout/rin)**(1./(nrad-1))
    for i in range(nrad):
        rad[i] = rin*dr**i

      
thet_min = 0.6981317008
thet_max = 1.57079637051
dthet  = (thet_max-thet_min)/(nthet-1)
thet = np.zeros(nthet)
for i in range(nthet):
    thet[i] = thet_min+i*dthet


phi_min = 0.0
phi_max = 6.28318548203
dphi  = (phi_max-phi_min)/(nphi-1)
phi = np.zeros(nphi)
for i in range(nphi):
    phi[i] = phi_min+i*dphi

if smallorbig == "small":
    opacity = np.genfromtxt('dustkappa_dust.inp.temp', skip_header=2, dtype=float)
elif smallorbig == "large":
    opacity = np.genfromtxt('dustkappa_dust.inp.imag', skip_header=2, dtype=float)
else:
    print 'ERROR in reading opacity file'
print " dustkappa read done."

opac_n, opac_m = opacity.shape
print opac_n, opac_m
for i in range(opac_n):
    if (opacity[i,0] < wavelength) and (opacity[i+1,0] > wavelength):
        print i
        interp_opac = i
        break
print interp_opac

opac_wave = opacity[interp_opac,1] + (opacity[interp_opac+1,1]-opacity[interp_opac,1])*(wavelength-opacity[interp_opac,0])/(opacity[interp_opac+1,0]-opacity[interp_opac,0])
print opac_wave




filename = "%s%s%s" % ("dust_density_", smallorbig, ".binp")
print filename

head = np.fromfile(filename, count=4, dtype=np.int64)
print "head[0]", head[0]
print "head[1]", head[1]
print "head[2]", head[2]
print "head[3]", head[3]

if head[2]!=(nrad*nthet*nphi):
    print ' ERROR'
    print ' Number of grid points in '+filename+' is different from that in amr_grid.inp'
    quit()

print " Reading the dust_temp file"

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

print " dust_density.inp read done."

density_3d = np.array(density_3d).reshape(nphi,nthet,nrad)
print " density_3d transformed to nphi,nthet,nrad"



# Calc the tau=1 surface along rad,phi (spherical coordinate)

tau_1 = np.zeros((nphi,nrad), dtype=float)
tau_all = np.zeros((nphi,nrad), dtype=float)
tau_1h = np.zeros((nphi,nrad), dtype=float)


for i in range(nphi):
    for j in range(nrad):
        tau = 0.0
        for k in range(nthet):
            tau = tau + opac_wave*density_3d[i,k,j]*rad[j]*np.sin(dthet)*1.496e13
            if (tau > 1.0):
                tau_1[i,j] = float(k)/float(nthet)
                r_v = rad[j]*math.cos(thet[k])
                r_r = rad[j]*math.sin(thet[k])
                h_r = 1.0 *(r_r/20.0**1.25)
                tau_h_tmp = r_v/h_r
                if tau_h_tmp > 100:
                    tau_1h[i,j] = 100
                elif tau_h_tmp < 0.1:
                    tau_1h[i,j] = 0.1
                else:
                    tau_1h[i,j] = tau_h_tmp


                #tau_h[i,j] = r_v/h_r
                #print "i,j,tau_1,k", i, j,k,tau_1[i,j]
                break
            elif (k == nthet-1) and (tau < 1):
                #print "optical thin"
                tau_1[i,j] = 1.0
                tau_h = 0.0

for i in range(nphi):
    for j in range(nrad):
        tau = 0.0
        for k in range(nthet):
            tau = tau + opac_wave*density_3d[i,k,j]*rad[j]*np.sin(dthet)*1.496e13
        if tau < 0.1:
            tau = 0.1
        if tau > 100:
            tau = 100
        tau_all[i,j] = tau

plt.contourf(rad, phi, tau_1, polar=True, cmap='jet')
#plt.contourf(rad, phi, tau_1, 15, alpha=.75, polar=True, cmap='jet')
#plt.contourf(rad, phi, surf_dens2d, 15, alpha=.75, polar=True, cmap='rainbow')
plt.xlabel('Radial Distance [AU]')
plt.ylabel('Azimuthal Angle')
cbar = plt.colorbar()
cbar.set_label('Nthet(tau=1) / Nthet(total) at %s micron' %wavelength, rotation=90)
plt.savefig('tau_1.png')

plt.clf()

plt.contourf(rad, phi, tau_all, norm=LogNorm(), cmap='jet')
plt.xlabel('Radial Distance [AU]')
plt.ylabel('Azimuthal Angle')
#plt.xlim([30,40])
cbar = plt.colorbar()
cbar.set_label('Tau to mid-plane at %s micron' %wavelength, rotation=90)
plt.savefig('tau_all.png')

plt.clf()

plt.contourf(rad, phi, tau_1h, norm=LogNorm(), cmap='jet')
plt.xlabel('Radial Distance [AU]')
plt.ylabel('Azimuthal Angle')
#plt.xlim([30,40])
cbar = plt.colorbar()
cbar.set_label('Tau to mid-plane at %s micron' %wavelength, rotation=90)
plt.savefig('tau_h.png')

#np.savetxt('au_phi_tau1_tauall.out', np.transpose([rad,oneradial]))


r_phi_t1_ta_th = np.zeros((nphi*nrad,5), dtype=float)
for i in range(nphi):
    for j in range(nrad):
        r_phi_t1_ta_th[i*nrad+j,0] = rad[j]
        r_phi_t1_ta_th[i*nrad+j,1] = phi[i]
        r_phi_t1_ta_th[i*nrad+j,2] = tau_1[i,j]
        r_phi_t1_ta_th[i*nrad+j,3] = tau_all[i,j]
        r_phi_t1_ta_th[i*nrad+j,4] = tau_1h[i,j]

np.savetxt('r_phi_t1_ta_th.out', np.transpose([r_phi_t1_ta_th[:,0],r_phi_t1_ta_th[:,1],r_phi_t1_ta_th[:,2],r_phi_t1_ta_th[:,3],r_phi_t1_ta_th[:,4]]))






