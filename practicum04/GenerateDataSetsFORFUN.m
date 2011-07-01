clear

% parameters
t0=0;
t1=10*365*24*3600;            % 10 years in sec

dur=t1-t0;                    % duration of the signal, spanning 10 years, in sec
Np=500;                       % number of data points in time ()
ut=t0+(0:Np-1)/Np*dur;        % evenly sampled data points
dt=dur/Np;                    % time resolution
f=0:1/dur:1/(2*dt);           % frequencies from DC to Nyquist
Nf=size(f,2);                 % Number of frequency bins


%create white noise data set  
for ll=1:2
  noisedata(ll,:)=randn(Np,1);    
end

% write white noise data set
for ll=1:2
    datafile=strcat('dataset1-',num2str(ll),'.txt');
    op=noisedata(ll,:);
    save(datafile,'op','-ASCII');
end

figure(1);plot(ut,noisedata)

%%%%%%%%%%%%% Add sinusoids to white noise

%create intrinsic white noise residual data for each pulsar    
for ll=1:2
  noisedata(ll,:)=noisedata(ll,:)+1/3*sin(2*pi*f(3)*ut)+sin(2*pi*f(100)*ut);
end

% write white noise data set
for ll=1:2
    datafile=strcat('dataset2-',num2str(ll),'.txt');
    op=noisedata(ll,:);
    save(datafile,'op','-ASCII');
end

figure(2);plot(ut,noisedata)

%%%%%%%%%%%%% Random data with noise burst

noiseburst=1.5*randn(Np,1); %burst of noise 1/5 of the duration

%create white noise data set  
for ll=1:2
  noisedata(ll,:)=randn(Np,1);    
end

for jj=300:400
    noisedata(1,jj)=(noisedata(1,jj)+noiseburst(jj));
    noisedata(2,jj)=(noisedata(2,jj)+noiseburst(jj));
end

% write white noise + burst data set
for ll=1:2
    datafile=strcat('dataset3-',num2str(ll),'.txt');
    op=noisedata(ll,:);
    save(datafile,'op','-ASCII');
end

figure(3);plot(ut,noisedata)





