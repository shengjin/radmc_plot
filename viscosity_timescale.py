#!/usr/bin/env python2
import math

import numpy as np

AU    = 1.4959787e13   # Astronomical unit [cm]
GG    = 6.67259e-8     # Gravitational constant
Msun  = 1.9891e33      # Solar mass [g]
Mstar = 1.00           # Stellar mass in Msun
alpha = 0.001

year  = 31556952.0    
h2r   = 0.10

r_chara = 10

Omega = math.sqrt(GG*Mstar*Msun/(r_chara*AU)**3.0)

vis   = alpha*h2r**2.0*Omega*(r_chara*AU)**2.0
tau_vis = (r_chara*AU)**2.0/vis/year

print "vis,tau_vis(year):"
print vis,tau_vis
