#!/usr/bin/env python2
import sys
import os
import linecache
import numpy as np
import math
import struct

#sys.path.append('/turquoise/users/sjin/python/')

import matplotlib as mpl
mpl.use('Agg')   # to avoid the window


rad_min = 0.24
rad_max = 15.6

scaling_radius = 10.0
#  // radius scale: LIN or LOG
scale  = "LIN" 


phi_min = 0.0
phi_max = 6.28318548203


filename1 = "%s%s" % ("./","/dust_density.cm.binp")
filename2 = "%s%s" % ("./","/dust_density.mm.binp")
amrname = "%s%s" % ("./", "/amr_grid.inp")
grid = linecache.getline(amrname, 6)
print filename1
print filename2
print amrname
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


'''
########################
###### Set up the grids

print "setting radial grid"

#  // define the radial grid
rin  = rad_min  * scaling_radius
rout = rad_max  * scaling_radius
rad = np.zeros(nrad)
if (scale == "LIN"):
    dr = (rout-rin)/(nrad-1) 
    for i in range(nrad):
        rad[i] = rin  + i*dr
if (scale == "LOG"):
    dr = (rout/rin)**(1./(nrad-1))
    for i in range(nrad):
        rad[i] = rin*dr**i

print "setting azimuth grid"
      

print "setting polar grid"

#  // define the azimuthal grid
dphi  = (phi_max-phi_min)/(nphi-1)
phi = np.zeros(nphi)
for i in range(nphi):
    phi[i] = phi_min+i*dphi

'''



###################
######## Start Function Declaration

######### ######### ######### ######### ######### ######### 
###### the function that read one temperature file (binary)

def readtemperature(filename):
    debug = False
 
    global temp_3d, nrad, nthet, nphi
    #temp_3d = np.array(temp_3d).reshape(nphi*nthet*nrad)

    head = np.fromfile(filename, count=4, dtype=np.int64)
    if debug:
	print "head[0]", head[0]
        print "head[1]", head[1]
        print "head[2]", head[2]
        print "head[3]", head[3]

    if head[2]!=(nrad*nthet*nphi):
        print ' ERROR'
        print ' Number of grid points in '+filename+' is different from that in amr_grid.inp'
        quit()

    if head[1] == 8:
        f = open(filename, "rb")  # reopen the file
        f.seek(32, os.SEEK_SET)  # seek
        temp_3dread = np.fromfile(f, count=-1, dtype=np.float64)
    elif head[1]==4:
        f = open(filename, "rb")  # reopen the file
        f.seek(32, os.SEEK_SET)  # seek
        temp_3dread = np.fromfile(f, count=-1, dtype=np.float32)
    else:
        print 'ERROR'
        print 'Unknown datatype in '+filename
    f.close()

    print " dust_temp reading: done."

    temp_3dread = np.array(temp_3dread).reshape(nphi,nthet,nrad)
    temp_3d = temp_3dread
    print " temp_3d transformed to nphi,nthet,nrad"

    return temp_3d



######## End Function Declaration
###################



################### the main program
################### ################### 

######  Declare the variables to hold the temperature data 
density_3d_1 = np.zeros((nphi,nthet,nrad), dtype=float)
density_3d_2 = np.zeros((nphi,nthet,nrad), dtype=float)


## Note: we use the upper nthet-nmodify layers of the first temperature data
readtemperature(filename1)
density_3d_1 = temp_3d
readtemperature(filename2)
density_3d_2 = temp_3d


########### write the new dust_density.binp
 
print "Write the new density data"


density_first = np.array(density_3d_1).reshape(nphi*nthet*nrad)
density_second = np.array(density_3d_2).reshape(nphi*nthet*nrad)
#np.savetxt('temp.dat', temp_3d) for ascii file
print "tranformed to 1 row array for write"

with open('dust_density.binp', 'wb') as f:
    f.write(struct.pack(4*'q',1,4,nphi*nthet*nrad,2))
    for num1 in density_first:
        f.write(struct.pack('f', num1))
    for num2 in density_second:
        f.write(struct.pack('f', num2))

