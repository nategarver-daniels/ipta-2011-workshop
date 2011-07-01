clear

T=365*24*3600; %one year

% physical and other constants
H0=2.4*10^(-18); % Hubble parameter in s^-1
c=3e8;           % speed of light in m/s
G=6.67e-11;      % Newton's constant in SI
Msun=2e30;       % Mass of teh sun in kg

% We would give them these cosmological functions valid for low redshifts
dz=0.01;      %redshift resolution for integrals
z=0:dz:2;     %redshift array that we integrate over
rz = c*H0^(-1)*z./(1+z/3.5042);    %proper distance in meters as function of redshift

%%%%%%% SMBBH %%%%%%%

% Here I use Alberto Sesana's by eye fit of the rate per unit redshift per 
% logarithmic interval of mass for the rate of SMBBH mergers
% (maybe we'd give this to Alberto Vecchio for his lecture)

Mstar=10^7.2 * Msun;  % parameter of Alberto's by eye fit
dlogM=0.1;            % logarithmic interval of mass for integration
logM=7:dlogM:100;     % mass integration range 
M=10.^(logM);
Cinv=sum( (M/Mstar).^0.7 .* exp(-(M/Mstar).^0.5 ))*dlogM; % rate normalization constant

% frequencies over which to compute the spectrum
logf=-10:0.1:-7;
f=10.^(logf);


% loop over frequencies (value of the spectrum for each frequency)
for ii=1:size(f,2)
    ii
    
    integrandINSP=0;
    for jj=2:size(z,2) % loop over redshifts
        
          for kk=1:size(M,2) % loop over masses

              % rate per unit redshift per logarithmic interval of
              % mass
              dRdzdlogM=3/4*1e-2 / Cinv * z(jj)^3 *(M(kk)/Mstar)^0.7 .* exp( -(M(kk)/Mstar)^0.5 )/(365*24*3600);
              
              Mt=2*M(kk);
              mu=M(kk)/2;
              
              % strain produced at frequency f by a source of mass M merging at redshift z, 
              hINSP=(G/c^2)^(5/6)*c^(1/6) *sqrt(pi/12) * (mu/rz(jj)) * Mt^(3/2)/mu^(1/2)* (pi*Mt*f(ii))^(-7/6) *(1+z(jj))^(-1/6); %from Eq. 44 of Thorne
                   
              integrandINSP=integrandINSP+dz*dlogM*hINSP^2*dRdzdlogM;   

          end 
    end
    
    % background spectrum in terms of Omega
    OmegaINSP(ii)=4*pi^2/3/H0^2 * f(ii)^3 *integrandINSP;
    
end

% Turn Omega into characteristic strain
hcSMBBH=sqrt(3/2/pi^2*H0^2*f.^(-2).*OmegaINSP);

figure;  loglog(f,hcSMBBH,'r'); grid on; hold on

% load Paul Demorest's futurisitic PTA frequency vs characteritic strain
load pta_2020.dat
f2=pta_2020(:,1);
hcPTA=pta_2020(:,2);

loglog(f2,hcPTA,'k')






