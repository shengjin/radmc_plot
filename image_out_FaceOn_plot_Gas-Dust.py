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

sizepix=image.sizepix_x

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


#!/usr/bin/env python2


# inclination of the disk
#incli = 51.0
incli = 0.0
#incli = 46.72

# position angle
# use it to scale the center Radius VS T_B
#posang = 61.0
posang = 138.02

# rotation angle used in the center ellipse to
rotang = posang
#rotang = 90-posang

# radius_wanted
n_rad = 100
inneravoid = 3
gap_width = 0.9 # AU

n_image = 0


dpc = 140.0  # distance in pc
n_image = 0  # the n_th image in image.out

# LkCa 15
#gaus_ma = 6.323864890469E-05    # maj axis of guassian in deg
#gaus_mi = 4.604819748137E-05    # min
#gaus_pa = 2.383443450928E+01    # PA of maj axis of Guass measure East from North

# HL Tau
gaus_ma = 8.367259158856E-06                                                  
gaus_mi = 5.335457002123E-06                                                  
gaus_pa = -1.758292236328E+02 

gaus_ma_arcsec = gaus_ma*3600
gaus_mi_arcsec = gaus_mi*3600
if gaus_pa > 0:
    gaus_pa_ctclk  = 180.0-gaus_pa  # PA in imConv
else:
    gaus_pa_ctclk  = -gaus_pa  # PA in imConv


fwhm = [gaus_ma_arcsec, gaus_mi_arcsec]

#debug = True
debug = False

AU = 1.496e13


#########################
#### print the gaussian
print "gaus_pa_ctclk"
print gaus_pa_ctclk
print "gaus_ma_arcsec, gaus_mi_arcsec"
print gaus_ma_arcsec, gaus_mi_arcsec
gaus_ma_rad = gaus_ma*math.pi/180.0
gaus_mi_rad = gaus_mi*math.pi/180.0
gaus_ma_AU = math.sin(gaus_ma_rad)*dpc*pc/AU
gaus_mi_AU = math.sin(gaus_mi_rad)*dpc*pc/AU
print "gaus_ma_AU, gaus_mi_AU"
print gaus_ma_AU, gaus_mi_AU


ImagePlot = ImageJyppix_scaled-Image_dustJyppix_scaled
#ImagePlot = image.imageJyppix[:,:,n_image]/dpc/dpc
print "Total flux mJy (original):"
print np.sum(ImagePlot, axis=None)*1000.0, "\n"


##########################################
# set the grids in pixel
########################################## 
# NOTE: in this case, nx=ny=m=n
# assuming the same number of grids in x,y axes
nx = image.nx
NN = nx
pix_x = np.linspace(-NN/2,NN/2,nx, endpoint=False)
ny = image.ny
pix_y = np.linspace(-NN/2,NN/2,ny, endpoint=False)

# set the grids in AU
im_x_au = (pix_x+0.5)*image.sizepix_x/au
im_y_au = (pix_y+0.5)*image.sizepix_y/au


############################################ 
# write the flux in Jy/pixel for debug
############################################ 
if debug:
    image_output = np.zeros(nx*nx)
    for i in range(nx):
        for j in range(nx):
            image_output[nx*i+j]=ImagePlot[i,j]
    np.savetxt('ImageJyppix_scaled', np.transpose(image_output))


##########################################
# define azimuthal extract function
#  could be ellipse or circle
##########################################

def azimuthal_Jy_avg(image, gap_au, gap_width, e_h, e_k, incli, rotang, sizepix):
    # input: image, gap_au, gap_width, e_h, e_k, incli, rotang, sizepix, nx
    # output: fin_vis
    ##################################
    # screen out the point in the ring at [i_in, i_out]

    au = 1.496e13

    fin0_x = []
    fin0_y = []
    fin0_vis = []
    
    # parameters for the center ellipse
    gap_min = gap_au - gap_width*0.5
    gap_max = gap_au + gap_width*0.5
    
    # assuming sizepix_x = sizepix_y
    e_a_grid_max = gap_max*au/sizepix  # long semi-axis
    e_a_grid_min = gap_min*au/sizepix  # long semi-axis
    
    inclination = math.radians(incli)
    e_b_grid_max = e_a_grid_max * cos(inclination)     # short semi-axis
    e_b_grid_min = e_a_grid_min * cos(inclination)     # short semi-axis
    rotang = math.radians(rotang)
    m = image.shape[0]
    # convert integer to float in order to make sure we find
    #    the center of the image.
    # image.out: 1) 20 grids
    #               python array: 0, 1, ..., 19
    #               center is at point 10
    #            2) 19 grids
    #               python array: 0, 1, ..., 18
    #               center is at point 9.5
    i_2_au = sizepix/au
    for ii in range(m):
        i=float(ii)
        for jj in range(m):
            j=float(jj)
            if ( ((i-e_h)*cos(rotang)+(j-e_k)*sin(rotang))**2/e_a_grid_max**2 + ((i-e_h)*sin(rotang)-( j-e_k)*cos(rotang))**2/e_b_grid_max**2 <= 1.0**2) and ( ((i-e_h)*cos(rotang)+(j-e_k)*sin(rotang))**2/e_a_grid_min**2 + ((i-e_h)*sin(rotang)-(j-e_k)*cos(rotang))**2/e_b_grid_min**2 > 1.0**2) :
                fin0_x.append(i)
                fin0_y.append(j)
                fin0_vis.append(image[ii,jj])
                #print nhalfpix, i, j
    fin_x = np.asarray(fin0_x)
    fin_y = np.asarray(fin0_y)
    fin_vis = np.asarray(fin0_vis)
    n_fin = fin_x.shape[0]
    if n_fin > 0: 
        np.savetxt('x_y_vis', np.transpose([fin_x,fin_y,fin_vis]))
    #total=np.sum(fin_vis), "\n"
    #avg = total/n_fin
    return fin_vis

##########################################
##########################################


nfloat = float(nx)
nhalfpix = nfloat/2
e_h = nhalfpix
e_k = nhalfpix

r_Jy = np.zeros([n_rad,2], dtype=float64)

for i in range(n_rad):
    gap_au = float(i)+inneravoid
    r_Jy[i,0] = gap_au
    avg = azimuthal_Jy_avg(ImagePlot, gap_au, gap_width, e_h, e_k, incli, rotang, sizepix)
    n_num = avg.shape[0]
    if n_num == 0:
        r_Jy[i,1] = 0
        print "WARNNING: no points found around ", r_Jy[i,0], " AU!!"
    else:
        f_n_num = float(n_num)
        #print avg
        total = np.sum(avg)
        #print total
        r_Jy[i,1] = total/f_n_num
        print r_Jy[i,0], r_Jy[i,1]
    
    

plt.xlabel('AU')
plt.ylabel("m Janskey / beam")
#plt.ylim(-180,180)
#plt.xlim(0,n_rad+inneravoid)
plt.plot(r_Jy[:,0], r_Jy[:,1]*1000)
#plt.pcolormesh(im_x_au, im_y_au, ImageJyppix_scaled*1000, vmin=1e-6, vmax=0.002,norm=LogNorm(),cmap='RdBu')
plt.savefig('azimuthalavg_fluxJy.png')
plt.clf()



#######################################################
#####  Make some Plot 
#########################################

plt.xlabel('AU')
plt.ylabel('AU')
plt.ylim(-180,180)
plt.xlim(-180,180)
plt.pcolormesh(im_x_au, im_y_au, ImagePlot*1000,cmap='RdBu')
#plt.pcolormesh(im_x_au, im_y_au, ImagePlot*1000,cmap='RdBu')
#plt.pcolormesh(im_x_au, im_y_au, ImageJyppix_scaled*1000, vmin=1e-6, vmax=0.002,norm=LogNorm(),cmap='RdBu')
cbar1=plt.colorbar()
cbar1.set_label("m Janskey / beam")
plt.title("Flux density")
plt.savefig('flux_density.png')
plt.clf()







