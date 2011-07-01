import numpy as np
import enthought.mayavi.mlab as mlab
#Some usefull websites:
# http://mathesaurus.sourceforge.net/matlab-numpy.html
# http://www.scipy.org/Tentative_NumPy_Tutorial
# https://www.cfa.harvard.edu/~jbattat/computer/python/science/idl-numpy.html
# http://code.enthought.com/projects/mayavi/docs/development/html/mayavi/mlab.html
#
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
mlab.show()
