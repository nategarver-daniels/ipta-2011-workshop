clear

T=365*24*3600; %one year in seconds

% physical and other constants
H0=2.4*10^(-18); % Hubble parameter in s^-1
c=3e8;           % speed of light in m/s
G=6.67e-11;      % Newton's constant in SI
Msun=2e30;       % Mass of the sun in kg

% We would give them these cosmological functions valid for low redshifts
dz=0.01;      %redshift resolution for integrals
z=0:dz:2;     %redshift array that we integrate over
rz = c*H0^(-1)*z./(1+z/3.5042);    %proper distance in meters as function of redshift