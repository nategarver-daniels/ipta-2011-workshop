#Practicum4B for students
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

directoryname='data/';

Np = 500             # number of data points
dur = 10 * 365 * 24 * 3600 # ten years of data
dt = dur/Np          # sampling time
f = np.arange(0,1/(2*dt),1/dur) # frequencies in steps of 1/dur up to Nyquist

####### load the first data set
noisedata1=np.empty([2,Np])

for ll in range(1,3):
    noisedata1[ll-1] = np.genfromtxt(directoryname+"dataset1-"+str(ll)+".txt")

meandata1 = np.mean(noisedata1[0])
std1 = np.std(noisedata1[0])
meandata2 = np.mean(noisedata1[1])
std2 = np.std(noisedata1[1])

correlation1 = 1/Np * np.sum((noisedata1[0] - meandata1) * (noisedata1[1] - meandata2))/(std1 * std2)

plt.figure(1)
plt.plot(range(Np),noisedata1[0])
plt.plot(range(Np),noisedata1[1])

noisedata2=np.empty([2,Np])

for ll in range(1,3):
    noisedata2[ll-1] = np.genfromtxt(directoryname+"dataset2-"+str(ll)+".txt")

# Find correlation for second dataset
meandata1 = np.mean(noisedata2[0])
std1 = np.std(noisedata2[0])
meandata2 = np.mean(noisedata2[1])
std2 = np.std(noisedata2[1])

correlation2 = 1/Np * np.sum((noisedata2[0] - meandata1) * (noisedata2[1] - meandata2))/(std1 * std2)

plt.figure(2)
plt.plot(range(Np),noisedata2[0],'r')
plt.plot(range(Np),noisedata2[1],'b')

# check for sinusoids in the data
y1 = np.fft.fft(noisedata2[0]) * dt 
y2 = np.fft.fft(noisedata2[1]) * dt 

plt.figure(3)
plt.plot(f,np.abs(y1[0:np.size(f)]),'r')
plt.plot(f,np.abs(y2[0:np.size(f)]),'b')

# add the data to see if sinusoids are correlated
# (i.e. in phase)
totaldata = noisedata2[0] + noisedata2[1]
y = np.fft.fft(totaldata) * dt

# Plot everything together
plt.figure(4)
plt.plot(f,np.abs(y1[0:np.size(f)]),'r')
plt.plot(f,np.abs(y2[0:np.size(f)]),'b')
plt.plot(f,np.abs(y[0:np.size(f)]),'g')

# Check that adding the data together in the frequency
# domain is equivalent
plt.figure(5)
plt.plot(f,np.abs(y1[0:np.size(f)]),'r')
plt.plot(f,np.abs(y2[0:np.size(f)]),'b')
plt.plot(f,np.abs(y1[0:np.size(f)] + y2[0:np.size(f)]),'g')

#Load the third data set
noisedata3=np.empty([2,Np])

for ll in range(1,3):
    noisedata3[ll-1] = np.genfromtxt(directoryname+"dataset3-"
+str(ll)+".txt")

meandata1 = np.mean(noisedata3[0])
std1 = np.std(noisedata3[0])
meandata2 = np.mean(noisedata3[1])
std2 = np.std(noisedata3[1])

correlation3 = 1/Np * np.sum((noisedata3[0] - meandata1) * (noisedata2[1] - meandata2))/(std1 * std2)

plt.figure(6)
plt.plot(range(Np),noisedata3[0],'r')
plt.plot(range(Np),noisedata3[1],'b')

# Try correlatin smaller chunks of the data

Nchunk = 25
correlation = np.zeros(Np - Nchunk)
for jj in range(Np - Nchunk):
    correlation[jj] = (1 / Nchunk) * np.sum(noisedata3[0,jj:Nchunk + jj]*noisedata3[1,jj:Nchunk + jj])

plt.figure(7)
plt.plot(range(Np - Nchunk),correlation*np.sqrt(Nchunk))

Nchunk = 50 
correlation = np.zeros(Np - Nchunk)
for jj in range(Np - Nchunk):
    correlation[jj] = (1 / Nchunk) * np.sum(noisedata3[0,jj:Nchunk + jj]*noisedata3[1,jj:Nchunk + jj])

plt.figure(8)
plt.plot(range(Np - Nchunk),correlation*np.sqrt(Nchunk))

Nchunk = 100 
correlation = np.zeros(Np - Nchunk)
for jj in range(Np - Nchunk):
    correlation[jj] = (1 / Nchunk) * np.sum(noisedata3[0,jj:Nchunk + jj]*noisedata3[1,jj:Nchunk + jj])

plt.figure(9)
plt.plot(range(Np - Nchunk),correlation*np.sqrt(Nchunk))

Nchunk = 150 
correlation = np.zeros(Np - Nchunk)
for jj in range(Np - Nchunk):
    correlation[jj] = (1 / Nchunk) * np.sum(noisedata3[0,jj:Nchunk + jj]*noisedata3[1,jj:Nchunk + jj])

plt.figure(10)
plt.plot(range(Np - Nchunk),correlation*np.sqrt(Nchunk))

Nchunk = 200 
correlation = np.zeros(Np - Nchunk)
for jj in range(Np - Nchunk):
    correlation[jj] = (1 / Nchunk) * np.sum(noisedata3[0,jj:Nchunk + jj]*noisedata3[1,jj:Nchunk + jj])

plt.figure(11)
plt.plot(range(Np - Nchunk),correlation*np.sqrt(Nchunk))

Nchunk = 250 
correlation = np.zeros(Np - Nchunk)
for jj in range(Np - Nchunk):
    correlation[jj] = (1 / Nchunk) * np.sum(noisedata3[0,jj:Nchunk + jj]*noisedata3[1,jj:Nchunk + jj])

plt.figure(12)
plt.plot(range(Np - Nchunk),correlation*np.sqrt(Nchunk))

plt.show()

