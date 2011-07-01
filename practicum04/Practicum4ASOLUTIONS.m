% script that generates bandpass FIR digital filter
clear

SamplingRate=64;           % 64Hz sampling rate
dur=20;                    % 20s duration filter

f=0:1/dur:SamplingRate/2;  %frequency vector from 0Hz to Nyquist in steps of 1/duration
Nf=length(f);              %length of frequency array

dt=1/SamplingRate;         % time between samples
t=0:dt:dur-dt;             % vector of times

% Frequency domain filter: Set frequency response to zero for all frequencies (for now)
f_domain_filter=zeros(length(f),1);

% low and high frequency cutoffs of the bandpass filter
f_low = 5;  % 5 Hz
f_high= 10; % 10 Hz

[i,j]=find((f < f_high) & (f > f_low));  % find frequencies lower than high freq. cutoff and higher than the low frequency cutoff
f_domain_filter(j)=1;      % set filter response to 1 at those frequencies

% plot filter in frequency domain for good measure
figure; plot(f,abs(f_domain_filter))

%this is a little bit of FFT magic: need to pack the FT properly so that ifft
%works correctly
for kk=1:Nf-2
      f_domain_filter(Nf+kk)=conj(f_domain_filter(Nf-kk));
end

% Inverse FT to generate FIR filter in time domain
t_domain_filter = ifft(f_domain_filter)/dt; %ifft includes (in matlab and I think python also) a factor of 1/N, so divide by dt to effectively multiply by df=1/T
    
% When generated this way the filter is symmetric around the beginning and
% end. This is due to the assumption of periodicity in the finite FFT. The end
% of the time series is the same as negative t. So
% we center the filter around zero.

% make filter symmetric around the center of the data stretch
len=length(t_domain_filter);
temp_filter(1:len/2)=t_domain_filter(len/2+1:len);
temp_filter(len/2+1:len)=t_domain_filter(1:len/2);

t_domain_filter=temp_filter;

t=-10:dt:10-dt;

figure; plot(t,t_domain_filter)

% If you zoom in on the filter you'll notice a lot of substructure. This is
% an example of the Gibbs phenomenon, where the sharp cutoff in frequency 
% domain filter leads to persistent ringing. A way to mitigate this effect is to add a
% window in the frequency domain to smoothly ramp up and down the band-pass
% filter

% find indices j consistent with low and high frequency cutoffs
[i,j]=find((f < f_high) & (f > f_low) );

% Set frequency response to zero for all frequencies
f_domain_filter=zeros(length(f),1);

% Hann window definition
hannwin = 1/2 * (1-cos(2*pi*(0:length(j)-1)/(length(j)-1)));

% Instead of a rectangular window like we used before use a Hann window
f_domain_filter(j)=hannwin;

% plot filter in frequency domain for good measure
figure; plot(f,abs(f_domain_filter))

% Now we regenerate the filter in the time domain (copy and paste from above)

% first pack the FT properly so that ifft works correctly
for kk=1:Nf-2
      f_domain_filter(Nf+kk)=conj(f_domain_filter(Nf-kk));
end

% then takle inverse FT to generate FIR filter in time domain
t_domain_filter = ifft(f_domain_filter)/dt; %ifft includes (in matlab and I think python also) a factor of 1/N, so divide by dt to effectively multiply by df=1/T

% finally make filter symmetric around the center of the data stretch
len=length(t_domain_filter);
temp_filter(1:len/2)=t_domain_filter(len/2+1:len);
temp_filter(len/2+1:len)=t_domain_filter(1:len/2);

t_domain_filter=temp_filter;

t=-10:dt:10-dt;

figure; plot(t,t_domain_filter)












