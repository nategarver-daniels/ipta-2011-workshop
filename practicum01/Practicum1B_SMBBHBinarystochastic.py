#Skeleton python script
from __future__ import division
import numpy as np
T = 365*24*3600 #one year in seconds

#Physical and other constants

H0 = 2.4E-18	#Hubble parameter in s^-1
c = 3E8 	#speed of light in m/s
G = 6.67E011	#Newton's constant in SI
Msun = 2E30	#Mass of the sun in kg

# Here are cosmological functions valid for _low_ redshifts

dz = 0.01	#redshift resolution for integrals
z=np.arange(0.0,2.0+dz,dz)
rz = c*H0**(-1)*z/(1 + z/3.5042) #proper distance in meters as function
				 #of redshift

print 3/4
