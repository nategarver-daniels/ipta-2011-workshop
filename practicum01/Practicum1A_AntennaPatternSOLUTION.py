import numpy as np
import enthought.mayavi.mlab as mlab
#Script that plots antenna pattern response of pulsar-earth system and LIGO

phi = np.linspace(0,2*np.pi,100)
theta = np.linspace(0,np.pi,100)

# Set up our grid

[phi,theta] = np.meshgrid(phi,theta)

# Radius of displayed sphere

r = 5

x = r*np.sin(theta)*np.cos(phi)
y = r*np.sin(theta)*np.sin(phi)
z = r*np.cos(theta)

#Plot Sphere

mlab.figure("Figure 1")
figure1 = mlab.mesh(x,y,z,representation='wireframe')
mlab.axes()
mlab.title("Figure 1")
mlab.outline()

#Make the radius a function of theta and phi
#antenna pattern response to a +/- polarized wave
#See Eq. 7 of Anholm et al. (Phys. Rev. D 79, 084030 (2009))

rp = np.absolute(.5*np.sin(theta)**2 * (np.cos(phi)**2 - 
	np.sin(phi)**2) / (1 + np.cos(theta)))

rp_mask = np.ma.masked_invalid(rp)

x = rp*np.sin(theta)*np.cos(phi)
x[np.isinf(x)] = np.nan
y = rp*np.sin(theta)*np.sin(phi)
y[np.isinf(y)] = np.nan
z = rp*np.cos(theta)
z[np.isinf(z)] = np.nan


mlab.figure("Figure 2")
figure2 = mlab.mesh(x,y,z,mask=np.ma.getmask(rp_mask),representation='wireframe')
mlab.axes()
mlab.title("Figure 2")
mlab.outline()

# Antenna pattern including frequency dependence (i.e. the pulsar term)
# (Eqs. 16 and 17 or Anholm et al.)

fL=10 #This is typical for pulsars
#Students should try this and notice a frequency dependence in  pictu

rp = np.absolute((np.exp(-2j*np.pi*fL*(1+np.cos(theta))) - 1)*
	.5*np.sin(theta)**2*(np.cos(phi)**2 - np.sin(phi)**2) 
	/ (1 + np.cos(theta)))

rp_mask = np.ma.masked_invalid(rp)

x = rp_mask*np.sin(theta)*np.cos(phi)
y = rp_mask*np.sin(theta)*np.sin(phi)
z = rp_mask*np.cos(theta)

mlab.figure("Figure 3")
figure3 = mlab.mesh(x,y,z,representation='wireframe')
mlab.axes()
mlab.title("Figure 3")
mlab.outline()
# EXCERCISE: At this point they do a Taylor expansion for theta=pi+delta, where delta
# is sufficiebtly small that the exponential can be expanded. This way they
# see where the frequency dependence in the response at the south pole is coming from 

# antenna pattern for LIGO (see Eq. 3 of Phys. Rev. D 79, 082002 (2009))
# this is kind of a cool paper, because it deals with non-Einstein tehories
# of gravity

# plus polarization with the polarization angle set to 0

psi = 0

rp = np.absolute(.5 * (1 + np.cos(theta)**2) * np.cos(2*phi) * np.cos(2*psi)
	- np.cos(theta)*np.sin(2*theta)*np.sin(2*psi))

rp_mask = np.ma.masked_invalid(rp)
x = rp_mask*np.sin(theta)*np.cos(phi)
y = rp_mask*np.sin(theta)*np.sin(phi)
z = rp_mask*np.cos(theta)

mlab.figure("Figure 4")
figure4 = mlab.mesh(x,y,z,representation='wireframe')
mlab.axes()
mlab.title("Figure 4")
mlab.outline()

mlab.show()
