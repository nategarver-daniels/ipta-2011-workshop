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




# read in first pulsar data set
residualdata1=np.empty([M,N])
for ll in range(M):
    residualdata1[ll] = np.genfromtxt(directoryname+"dataset1-"+str(ll+1)+".txt")
    # You may want to plot the data after you've read it in


## (B) Write some code to computes and stores the correlation for each pulsar pair.





## (C) Then plot the correlation for each pair from (B) versus the angular separation you
## calculated in (A)





## (D) Finally calculate Eq. (4) of Ref. 1, and the significance (multiplying by the
## square root of the number of pulsar pairs) using your results from (A) and (B)



# Once everything is working: copy, paste, and edit appropriately for data sets 2 and 3

