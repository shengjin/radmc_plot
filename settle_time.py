#!/usr/bin/python2 

import numpy as np
#sys.path.append('/turquoise/users/sjin/python/')

import matplotlib as mpl
mpl.use('Agg')   # to avoid the window

from matplotlib import colors
from pylab import setp
import matplotlib.pyplot as plt


font = {'family' : 'serif', #monospace
        'weight' : 'normal', #bold
        'size'   : 12,
        }

font1 = {'family' : 'monospace', #monospace
        'weight' : 'bold', #bold
        'size'   : 15,
        }

mpl.rc('font', **font)

filename_gas = "%s%s" % ("gas","_surface_density.ascii")

PI=3.14159265358979323846264338328
TWOPI=6.28318530717958647692528676656
PIhalf=1.57079632679489661923132169164

AU2CM  = 1.4959787E13

Mstar = 0.55*1.9891E33
G = 6.67259E-8

Sec2Yr=31557600.0

theta_r_rho_gas = np.genfromtxt(filename_gas, dtype=float)

r=theta_r_rho_gas[:,1]
theta=theta_r_rho_gas[:,0]
Sigma_gas=theta_r_rho_gas[:,2]

x_gas=r*np.cos(theta)/AU2CM
y_gas=r*np.sin(theta)/AU2CM
#stokesnumber = (1.0/rho)*0.3*1.57079637051
Omega_k = (G*Mstar/(r**3.0))**0.5
#Omega_k = (G*Mstar/(r**3.0)) #**0.5
v_k = (G*Mstar/r)**0.5
rho_d = 1.0 # 3.0 g cm^-3
#a = 0.0001 # diameter of 1 micron grains in cm g s^-1
a = 0.015 # diameter of 0.15 mm grains in cm g s^-1

h_hydro = 0.05*(r/10/AU2CM)**0.285*r
h_radmc = 1.0*(r/20/AU2CM)**1.25*AU2CM
#h_hydro = 0.05*(r/10/AU2CM)**0.285
#h_radmc = 1.0*(r/20/AU2CM)**1.25/r
h = h_hydro
#h = h_radmc
#z = 0.1*h # 
z = 1.0*h # 
rho_z = np.exp(-z**2.0/(2*h**2.0))*Sigma_gas/(h*TWOPI**0.5)

cs = h*Omega_k
#cs = h/r*v_k
v_bar = (4/PIhalf)**0.5*cs

v_settle = (Omega_k**2.0/v_bar)*rho_d*a*z/rho_z

t_settle1 = z/np.abs(v_settle)
t_settle1 = t_settle1/Sec2Yr
t_settle2 = (2/PI)*np.exp(-z**2.0/(2*h**2.0))*Sigma_gas/(Omega_k*rho_d*a)
t_settle2 = t_settle2/Sec2Yr

t_fric = rho_d*a/rho_z/v_bar

plt, axes = plt.subplots(nrows=1, ncols=2, figsize=(11, 6.5), dpi=180, sharex=True, sharey=True)

axisadd=10

im=axes[0].scatter(x_gas, y_gas, c=t_settle1, edgecolors='none', marker="o", s=(r/AU2CM/20), vmin=t_settle1.min(), vmax=t_settle1.max(), cmap="jet")
axes[0].axis([x_gas.min()-axisadd, x_gas.max()+axisadd, y_gas.min()-axisadd, y_gas.max()+axisadd])
axes[0].text(-158, 142, 'z = h_gas', fontdict=font1)
axes[0].set_xlabel(r'$\mathrm{Distance \hspace{0.4} [ \hspace{0.3} AU \hspace{0.3}]}$', rotation=0, fontsize=16,  labelpad=5)
axes[0].set_ylabel(r'$\mathrm{Distance \hspace{0.4} [ \hspace{0.3} AU \hspace{0.3}]}$', rotation=90, fontsize=16,  labelpad=0)
#fig.colorbar(cs, ax=ax, shrink=0.9)
#ticks_at1 = [0,20,40,60,80,100,120,140]
#cbar1=plt.colorbar(im, ax=axes[0], ticks=ticks_at1, orientation='horizontal', shrink=0.8)
cbar1=plt.colorbar(im, ax=axes[0], orientation='horizontal', format='%.0e', shrink=1.2)
cbar1.ax.tick_params(labelsize=10)
cbar1.set_label(r'$\mathrm{Settling \hspace{0.4} Time \hspace{0.4} [ \hspace{0.3} Year \hspace{0.3}]}$', rotation=0, fontsize=14,  labelpad=0)


im=axes[1].scatter(x_gas, y_gas, c=v_settle, edgecolors='none', marker="o", s=(theta_r_rho_gas[:,1]/AU2CM/30), vmin=v_settle.min(), vmax=v_settle.max(), cmap="jet")
axes[1].axis([x_gas.min()-axisadd, x_gas.max()+axisadd, y_gas.min()-axisadd, y_gas.max()+axisadd])
axes[1].text(-145, 142, 'rho_d = 1 g/cm^-3', fontdict=font1)
axes[1].set_xlabel(r'$\mathrm{Settling \hspace{0.4} Velocity [ \hspace{0.3} cm/s \hspace{0.3} ]}$', rotation=0, fontsize=16,  labelpad=5)
#ticks_at2 = [0,0.5,1.0,1.5,2.0,2.5]
#cbar2=plt.colorbar(im, ax=axes[1], ticks=ticks_at2, orientation='horizontal', format='%.1f', shrink=0.8)
cbar2=plt.colorbar(im, ax=axes[1], orientation='horizontal', format='%.1f', shrink=1.2)
cbar2.ax.tick_params(labelsize=12)
cbar2.set_label(r'$\mathrm{Settling \hspace{0.4} Velocity \hspace{0.4} [ \hspace{0.3}cm/s \hspace{0.3}]}$', rotation=0, fontsize=14)


#im=axes[2].scatter(x_gas, y_gas, c=h, edgecolors='none', marker="o", s=(theta_r_rho_gas[:,1]/AU2CM/30), vmin=h.min(), vmax=h.max(), cmap="jet")
#im=axes[2].scatter(x_gas, y_gas, c=v_bar, edgecolors='none', marker="o", s=(theta_r_rho_gas[:,1]/AU2CM/30), vmin=v_bar.min(), vmax=v_bar.max(), cmap="jet")
#im=axes[2].scatter(x_gas, y_gas, c=Omega_k, edgecolors='none', marker="o", s=(theta_r_rho_gas[:,1]/AU2CM/30), vmin=Omega_k.min(), vmax=Omega_k.max(), cmap="jet")
#im=axes[2].scatter(x_gas, y_gas, c=t_fric, edgecolors='none', marker="o", s=(theta_r_rho_gas[:,1]/AU2CM/30), vmin=t_fric.min(), vmax=t_fric.max(), cmap="jet")
#im=axes[2].scatter(x_gas, y_gas, c=rho_z, edgecolors='none', marker="o", s=(theta_r_rho_gas[:,1]/AU2CM/30), vmin=rho_z.min(), vmax=rho_z.max(), cmap="jet")
#axes[2].axis([x_gas.min()-axisadd, x_gas.max()+axisadd, y_gas.min()-axisadd, y_gas.max()+axisadd])
#axes[2].text(-145, 142, 't_fric', fontdict=font1)
#axes[2].set_xlabel(r'$\mathrm{Friction \hspace{0.4} Time [ \hspace{0.3} s \hspace{0.3} ]}$', rotation=0, fontsize=16,  labelpad=5)
##ticks_at2 = [0,0.5,1.0,1.5,2.0,2.5]
##cbar2=plt.colorbar(im, ax=axes[1], ticks=ticks_at2, orientation='horizontal', format='%.1f', shrink=0.8)
#cbar3=plt.colorbar(im, ax=axes[2], orientation='horizontal', format='%.0e', shrink=1.2)
#cbar3.ax.tick_params(labelsize=8)
#cbar3.set_label(r'$\mathrm{Settling \hspace{0.4} Velocity \hspace{0.4} [ \hspace{0.3}cm/s \hspace{0.3}]}$', rotation=0, fontsize=14)


plt.tight_layout()

plt.savefig('Tsett_Vsett.png')

quit()












