clear

directoryname='data/';

Np=500;             % number of data points
dur=10*365*24*3600; % ten years of data
dt=dur/Np;          % sampling time
f=0:1/dur:1/(2*dt); % frequencies in steps of 1/dur up to Nyquist

%%%%%%% load the first data set
for ll=1:2
    datafile=strcat(directoryname,'dataset1-',num2str(ll),'.txt');
    noisedata(ll,:)=load(datafile);
end



