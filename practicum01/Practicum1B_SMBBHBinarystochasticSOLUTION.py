#Skeleton python script
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
np.seterr(divide='ignore',over='ignore')
T = 365 * 24 * 3600 #one year in seconds

#Physical and other constants

H0 = 2.4 * 10**(-18)	#Hubble parameter in s^-1
c = 3E8 	#speed of light in m/s
G = 6.67E-11	#Newton's constant in SI
Msun = 2E30	#Mass of the sun in kg


# Here are cosmological functions valid for _low_ redshifts

dz = 0.01	#redshift resolution for integrals
z = np.arange(0.0,2.0 + dz,dz)
rz = c * H0**(-1) * z * (1 + z / 3.5042)**-1 #proper distance in meters as function
				 #of redshift
####### SMBBH #######

# Here I use Alberto Sesana's by eye fit of the rate per unit redshift per 
# logarithmic interval of mass for the rate of SMBBH mergers
# (maybe we'd give this to Alberto Vecchio for his lecture)

Mstar = 10**7.2 * Msun	#Parameter of Alberto's by eye fit
dlogM = 0.1
logM = np.arange(7.0,100 + dlogM,dlogM)
M = 10**logM
Cinv = ((M/Mstar)**0.7 * np.exp(-(M/Mstar)**0.5)).sum(axis=0) * dlogM #rate
#normalization constant

#Frequencies over which to computer the spectrum

logf = np.arange(-10,-6.9,.1)
f = 10**logf
#Loop over frequencies (value of the spectrum for each frequency)


OmegaINSP=np.zeros(f.shape[0])
for ii in range(0,np.size(f)):
    integrandINSP = 0.
    for jj in range(1,np.size(z)):
        #Rate per unit redshift per logarithmic interval of mass
        for kk in range(0,np.size(M)-1):
            dRdzdlogM = 3/4 * 1E-2 / Cinv * z[jj]**3 * (M[kk] / Mstar)**0.7 * np.exp(-(M[kk] / Mstar)**0.5) / (365 * 24 * 2600)
            Mt = 2 * M[kk]
            mu = M[kk] / 2
            # Strain produced at freq. f by a source of mass M merging at 
            # redshift z

            hINSP=(G / c**2)**(5/6)*c**(1/6)*np.sqrt(np.pi/12)*(mu / rz[jj]) * Mt**(3/2)/mu**(1/2)*(np.pi*Mt*f[ii])**(-7/6)*(1+z[jj])**(-1/6)
            integrandINSP = integrandINSP + dz * dlogM * hINSP**2 * dRdzdlogM
    OmegaINSP[ii] = 4*np.pi**2 / 3 / H0**2 * f[ii]**3 * integrandINSP
hcSMBBH=np.sqrt(3/2/np.pi**2*H0**2*f**-2*OmegaINSP)

#Load data from pta file
data = np.genfromtxt('pta_2020.dat',unpack=True,usecols=(0,1),skip_header=1)
plt.figure(1)
plt.loglog(f,hcSMBBH)
plt.loglog(data[0],data[1])
plt.grid(which='Minor')
plt.show()
