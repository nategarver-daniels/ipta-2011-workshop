# script that generates bandpass FIR digital filter
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


SamplingRate=64;           # 64Hz sampling rate
dur=20;                    # 20s duration filter

f=np.arange(0,SamplingRate/2 + 1/dur,1/dur);  #frequency vector from 0Hz to Nyquist in steps of 1/duration
Nf=np.size(f);              #length of frequency array

dt=1/SamplingRate;         # time between samples
t=np.arange(0,dur,dt);             # vector of times

# Frequency domain filter: Set frequency response to zero for all frequencies (for now)
f_domain_filter=np.zeros(Nf);

# low and high frequency cutoffs of the bandpass filter
f_low = 5;  # 5 Hz
f_high= 10; # 10 Hz

i = np.where((f < f_high) & ( f > f_low))
# find frequencies lower than high freq. cutoff and 
#higher than the low frequency cutoff

f_domain_filter[i]=1;      # set filter response to 1 at those frequencies

# plot filter in frequency domain for good measure
figure1 = plt.figure(1)
plt.plot(f,np.abs(f_domain_filter))
#this is a little bit of FFT magic: need to pack the FT properly so that ifft
#works correctly

#print np.size(f_domain_filter[-2:0:-1])
#f_domain_filter_conj = f_domain_filter[-2:0:-1].conj()
f_domain_filter=np.concatenate((f_domain_filter,f_domain_filter[-2:0:-1].conj()))

# Inverse FT to generate FIR filter in time domain
t_domain_filter = np.fft.ifft(f_domain_filter)/dt;
#ifft includes (in matlab and I think python also) 
#a factor of 1/N, so divide by dt to effectively multiply by df=1/T
# When generated this way the filter is symmetric around the beginning and
# end. This is due to the assumption of periodicity in the finite FFT. The end
# of the time series is the same as negative t. So
# we center the filter around zero.
#print t_domain_filter 
# make filter symmetric around the center of the data stretch
length = np.size(t_domain_filter);
temp_filter = np.concatenate((t_domain_filter[length/2::1],t_domain_filter[:length/2:1]))
t_domain_filter = temp_filter
###temp_filter(1:len/2)=t_domain_filter(len/2+1:len);
###temp_filter(len/2+1:len)=t_domain_filter(1:len/2);
###t_domain_filter=temp_filter;

t=np.arange(-10,10,dt);
# If you zoom in on the filter you'll notice a lot of substructure. This is
# an example of the Gibbs phenomenon, where the sharp cutoff in frequency 
# domain filter leads to persistent ringing. A way to mitigate this effect is to add a
# window in the frequency domain to smoothly ramp up and down the band-pass
# filter.  This is what you'll do for this part of the practicum.

figure2 = plt.figure(2)
plt.plot(t,t_domain_filter.real)

f_domain_filter = np.zeros(Nf)

hannwin = 1/2 * (1-np.cos(2*np.pi*np.arange(np.size(i) + 1)/(np.size(i))))

f_domain_filter[i]=hannwin

figure3 = plt.figure(3)
plt.plot(f,np.abs(f_domain_filter))

f_domain_filter=np.concatenate((f_domain_filter,f_domain_filter[-2:0:-1].conj()))

t_domain_filter = np.fft.ifft(f_domain_filter)/dt;

t_domain_filter=np.roll(t_domain_filter,(Nf - 1))

figure4 = plt.figure(4)
plt.plot(t,t_domain_filter.real)
plt.show()

