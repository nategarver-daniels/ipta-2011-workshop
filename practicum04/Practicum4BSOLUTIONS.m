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

meandata1=mean(noisedata(1,:));
std1=std(noisedata(1,:));
meandata2=mean(noisedata(2,:));
std2=std(noisedata(2,:));

% find correlation for first dataset
correlation1=1/Np*sum((noisedata(1,:)-meandata1).*(noisedata(2,:)-meandata2))/(std1*std2)     
figure; plot(1:Np,noisedata)

%%%%%% load the second data set
for ll=1:2
    datafile=strcat(directoryname,'dataset2-',num2str(ll),'.txt');
    noisedata(ll,:)=load(datafile);
end

% find correlation for second dataset
meandata1=mean(noisedata(1,:));
std1=std(noisedata(1,:));
meandata2=mean(noisedata(2,:));
std2=std(noisedata(2,:));

correlation2=1/Np*sum((noisedata(1,:)-meandata1).*(noisedata(2,:)-meandata2))/(std1*std2)     
figure; plot(1:Np,noisedata)

figure;plot(1:Np,noisedata)

% check for sinusoids in the data
y1=fft(noisedata(1,:))*dt;
y2=fft(noisedata(2,:))*dt;

figure;plot(f,abs(y1(1:length(f)))); hold on; plot(f,abs(y2(1:length(f))),'r')

% add the data to see if sinusoids are correlated (i.e. in phase)
totaldata=noisedata(1,:)+noisedata(2,:);
y=fft(totaldata)*dt;

% plot everything together
figure;plot(f,abs(y1(1:length(f)))); hold on; plot(f,abs(y2(1:length(f))),'r');
plot(f,abs(y(1:length(f))),'g')

% check that adding the data together in the frequency domain is equivalent
figure;plot(f,abs(y1(1:length(f)))); hold on; plot(f,abs(y2(1:length(f))),'r');
plot(f,abs(y1(1:length(f))+y2(1:length(f))),'g')

%%%%%% load third data set
for ll=1:2
    datafile=strcat(directoryname,'dataset3-',num2str(ll),'.txt');
    noisedata(ll,:)=load(datafile);
end
% find correlation for second dataset
meandata1=mean(noisedata(1,:));
std1=std(noisedata(1,:));
meandata2=mean(noisedata(2,:));
std2=std(noisedata(2,:));

correlation3=1/Np*sum((noisedata(1,:)-meandata1).*(noisedata(2,:)-meandata2))/(std1*std2)     
figure;plot(1:Np,noisedata)

% try correlating smaller chunks of the data
Nchunk=25;
for jj=1:Np-Nchunk
    correlation(jj)=1/Nchunk*sum(noisedata(1,jj:Nchunk+jj).*noisedata(2,jj:Nchunk+jj));
end
figure;plot(correlation*sqrt(Nchunk))
clear correlation


Nchunk=50;
for jj=1:Np-Nchunk
    correlation(jj)=1/Nchunk*sum(noisedata(1,jj:Nchunk+jj).*noisedata(2,jj:Nchunk+jj));
end
figure;plot(correlation*sqrt(Nchunk))
clear correlation

Nchunk=100;
for jj=1:Np-Nchunk
    correlation(jj)=1/Nchunk*sum(noisedata(1,jj:Nchunk+jj).*noisedata(2,jj:Nchunk+jj));
end
figure;plot(correlation*sqrt(Nchunk))
clear correlation

Nchunk=150;
for jj=1:Np-Nchunk
    correlation(jj)=1/Nchunk*sum(noisedata(1,jj:Nchunk+jj).*noisedata(2,jj:Nchunk+jj));
end
figure;plot(correlation*sqrt(Nchunk))
clear correlation


Nchunk=200;
for jj=1:Np-Nchunk
    correlation(jj)=1/Nchunk*sum(noisedata(1,jj:Nchunk+jj).*noisedata(2,jj:Nchunk+jj));
end
figure;plot(correlation*sqrt(Nchunk))
clear correlation

Nchunk=250;
for jj=1:Np-Nchunk
    correlation(jj)=1/Nchunk*sum(noisedata(1,jj:Nchunk+jj).*noisedata(2,jj:Nchunk+jj));
end
figure;plot(correlation*sqrt(Nchunk))
clear correlation


