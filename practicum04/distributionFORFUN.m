clear

dr=0.001;

r=-1:dr:1;

Np=500;

p=(1-r.^2).^((Np-2)/2); % un-normalised pdf

p=p/(sum(p)*dr); % normalised pdf

figure;plot(r,p)
