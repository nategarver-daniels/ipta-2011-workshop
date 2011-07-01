#Practicum5 starter script
from __future__ import division
import numpy as np 
import matplotlib.pyplot as plt

# directory with all data
directoryname='data/';

# load pulsar sky location file
phi,theta = np.loadtxt(directoryname+"SkyPositions.txt", unpack=True, usecols=[0,1])

# Number of pulsars
M=20

# Number of data points in each file
N=500

## (A) Write some code that loops over pulsar pairs and determines
## the angular separation for each pair, as well as the expected value of the HD curve
## for each pair
xi = np.zeros(N)
angularseparation = np.zeros(N)
## Here, I've over-provisioned xi and angularseparation, because it's speedier
## to fill and then truncate a zero-initialized array that to append to an 
## existing array in NumPy

jj = 0 
for ll in range(M):
    phati = np.array([np.cos(phi[ll])*np.sin(theta[ll]),np.sin(phi[ll])*np.sin(theta[ll]),np.cos(theta[ll])])
    for kk in range(ll+1,M):
        phatj = np.array([np.cos(phi[kk])*np.sin(theta[kk]),np.sin(phi[kk])*np.sin(theta[kk]),np.cos(theta[kk])])
        angularseparation[jj] = np.arccos(np.sum(phati * phatj))
        xip = (1 - np.sum(phati * phatj))/2
        xi[jj] = 3 * (1/3 + xip * (np.log(xip) - 1/6))/2
        jj = jj+1
# As mentioned earlier, now we truncate xi and angularseparation
xi = np.trim_zeros(xi)
angularseparation = np.trim_zeros(angularseparation)


# read in first pulsar data set
residualdata1=np.empty([M,N])
for ll in range(M):
    residualdata1[ll] = np.genfromtxt(directoryname+"dataset1-"+str(ll+1)+".txt")
    # You may want to plot the data after you've read it in


plt.figure(1)
for line in range(M):
    plt.plot(range(N),residualdata1[line])

## (B) Write some code to computes and stores the correlation for each pulsar pair.
## We know correlation will be the same size as angularseparation, so lets
## size based on that (no trimming)
correlation = np.zeros(np.size(angularseparation))

jj = 0
for ll in range(M):
    for kk in range(ll+1,M):
        correlation[jj] = 1 / N * np.sum(residualdata1[ll] * residualdata1[kk])
        jj = jj + 1

print np.shape(residualdata1)
print correlation
## (C) Then plot the correlation for each pair from (B) versus the angular separation you
## calculated in (A)
plt.figure(2)
plt.scatter(angularseparation,correlation)

## (D) Finally calculate Eq. (4) of Ref. 1, and the significance (multiplying by the
## square root of the number of pulsar pairs) using your results from (A) and (B)

Npairs = M * (M - 1)/2

# Calculate the coherence
#
# First calculate the mean
meanr = np.mean(correlation)
meanxi = np.mean(xi)

# Then calculate standard deviation

stdr = np.sqrt(1/Npairs * np.sum((correlation - meanr)**2))
stdxi = np.sqrt(1/Npairs * np.sum((xi - meanxi)**2))

rho = 1/Npairs * np.sum((correlation - meanr) * (xi - meanxi)) / (stdr * stdxi)

Significance = rho * np.sqrt(Npairs)

print "rho: ",rho
print "Significance: ",Significance

# Once everything is working: copy, paste, and edit appropriately for data sets 2 and 3

