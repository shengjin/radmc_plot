#!/usr/bin/env python2

import numpy as np

import matplotlib as mpl
mpl.use('Agg')  # to avoid window

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from readimage import *
from natconst import *
from matplotlib import cm
import math
from math import sin, cos


#debug = True
debug = False


# inclination of the disk
incli = 46.72
# position angle
# use it to scale the center Radius VS T_B
posang = 138.02

# rotation angle used in the center ellipse to
# manually increase the center flux
rotang = 90-138.02

###################### ###################### 
########### Read image and scale the flux
###################### ###################### 

image = readImage('image_hydroV.out')
dpc = 140.0  # distance in pc
n_image = 0  # the n_th image in image.out
ImageJyppix_scaled = image.imageJyppix[:,:,n_image]/dpc/dpc
print "Total flux mJy (original):"
print np.sum(ImageJyppix_scaled, axis=None)*1000.0, "\n"

image_dust = readImage('image_3361dust.out')
Image_dustJyppix_scaled = image_dust.imageJyppix[:,:,n_image]/dpc/dpc

##########################################
# set the grids in pixel
########################################## 

#################################
# NOTE: in this case, nx=ny=m=n
# assuming the same number of grids in x,y axes
nx = image.nx
NN = nx
pix_x = np.linspace(-NN/2,NN/2,nx, endpoint=False)
fft_dx = pix_x[1] - pix_x[0]
print pix_x
ny = image.ny
pix_y = np.linspace(-NN/2,NN/2,ny, endpoint=False)
fft_dy = pix_y[1] - pix_y[0]

n_mid=(NN-1)/2
x_to_au = (pix_x+0.5)*image.sizepix_x/au
#np.savetxt('T_B.out', np.transpose([T_B[n_mid,:], x_to_au]))

############################################ 
# write the flux in Jy/pixel for debug
############################################ 
if debug:
    image_output = np.zeros(nx*nx)
    for i in range(nx):
        for j in range(nx):
            image_output[n*i+j]=ImageJyppix_scaled[i,j]
    np.savetxt('ImageJyppix_scaled', np.transpose(image_output))



line1 = np.zeros(nx)
line2 = np.zeros(nx)
line3 = np.zeros(nx)
line4 = np.zeros(nx)

#np.savetxt('Fl.out', ImageJyppix_scaled[n_mid,:]/dpc/dpc)

line1 = ImageJyppix_scaled[n_mid,:]-Image_dustJyppix_scaled[n_mid,:]
line2 = ImageJyppix_scaled[n_mid+1,:]-Image_dustJyppix_scaled[n_mid+1,:]
line3 = ImageJyppix_scaled[:,n_mid]-Image_dustJyppix_scaled[:,n_mid]
line4 = ImageJyppix_scaled[:,n_mid+1]-Image_dustJyppix_scaled[:,n_mid+1]
np.savetxt('line1.out', np.transpose([line1, x_to_au]))
np.savetxt('line2.out', np.transpose([line2, x_to_au]))
np.savetxt('line3.out', np.transpose([line3, x_to_au]))
np.savetxt('line4.out', np.transpose([line4, x_to_au]))
 
############################################ 
############ info printting
############################################ 

print "pixel-size (x) in cm :"
print "     ", image.sizepix_x
print "pixel-size (y) in cm :"
print "     ", image.sizepix_y
print "distance from the source (cm), PC x dpc(nPCs) :"
print "     ", dpc*3.0856775e18
print "radian in x :"
print "     ", image.sizepix_x/dpc/3.0856775e18
print "radian in y :"
print "     ", image.sizepix_y/dpc/3.0856775e18, "\n"


plt.xlim(3,150)
plt.loglog(x_to_au,line1,label='flux_y1', color='black', lw=1, ls='-.')
plt.loglog(x_to_au,line2,label='flux_y2', color='black', lw=1, ls='-.')
plt.loglog(x_to_au,line3,label='flux_x1', color='black', lw=1, ls='-.')
plt.loglog(x_to_au,line4,label='flux_x2', color='black', lw=1, ls='-.')
plt.savefig('radial_flux.png')
plt.clf()

##########################
# make a plot
plt.pcolormesh(pix_x, pix_y, ImageJyppix_scaled-Image_dustJyppix_scaled, cmap=cm.gray)
plt.colorbar()
#plt.ylim(-180,180)
grids = "%s%s%s%s" % (nx+1, " x ", ny+1, " grids")
plt.title(grids)
plt.savefig('ori.png')
plt.clf()


#####################
# set the grids in AU
# change
im_x_au = (pix_x)*image.sizepix_x/au
im_y_au = (pix_y)*image.sizepix_y/au


#######################################################
#####  Make some Plot 
#########################################

plt.xlabel('AU')
plt.ylabel('AU')
plt.ylim(-100,100)
plt.xlim(-100,100)
plt.pcolormesh(im_x_au, im_y_au, ImageJyppix_scaled-Image_dustJyppix_scaled,cmap='RdBu')
cbar1=plt.colorbar()
cbar1.set_label("Janskey/pixel")
plt.title("Flux density")
plt.savefig('flux_density.png')
plt.clf()



