clear

%%% Random data set %%%
%%%%% Number of pulsars
%Npulsars = 20;
%%%%% Pulsar sky-locations
%phi=2*pi*rand(Npulsars,1);     %random ascension in radians (0-2pi)
%theta=pi*(1-rand(Npulsars,1)); %random declination in radians (0-pi)
% save the sky positions file
%for ll=1:Npulsars
%    op1(ll,1)=phi(ll);
%    op1(ll,2)=theta(ll);
%end
%skypositionsfile='SkyPositions.txt';
%save(skypositionsfile,'op1','-ASCII');

%%%% Fixed data set %%%%
load SkyPositions.txt
Npulsars=length(SkyPositions);
phi=SkyPositions(:,1);
theta=SkyPositions(:,2);

% parameters
t0=0;
t1=10*365*24*3600;            % 10 years in sec

dur=t1-t0;                    % duration of the signal, spanning 10 years, in sec
Np=500;                       % number of data points in time ()
ut=t0+(0:Np-1)/Np*dur;        % evenly sampled data points
dt=dur/Np;                    % time resolution
f=0:1/dur:1/(2*dt);           % frequencies from DC to Nyquist
Nf=size(f,2);                 % Number of frequency bins
ResidualRMS=1e-7;             % RMS of intrinsic white residuals in seconds

% Current Hubble scale in 1/seconds
H0=2.27e-18;

% Stochastic background parameters in terms of Omega
alpha=0;
Omega_alpha = 1e-9;  %Here I set the value of the background to unity because I want S
Omega=Omega_alpha*f.^alpha;

%create intrinsic white noise residual data for each pulsar    
for ll=1:Npulsars
  residualdata(ll,:)=ResidualRMS*randn(Np,1);    
end

% Create random frequency series from zero mean, unit variance, Gaussian
% distributions
for ll=1:Npulsars
    wlocal(ll,:)=randn(Nf,1)+i*randn(Nf,1);      
end

% create overlap reduction function matrix (Npulsars x Npulsars)
for ll=1:Npulsars
    
      phati(1)=cos(phi(ll))*sin(theta(ll));
      phati(2)=sin(phi(ll))*sin(theta(ll));
      phati(3)=cos(theta(ll));
    
    for kk=1:Npulsars
      
      phatj(1)=cos(phi(kk))*sin(theta(kk));
      phatj(2)=sin(phi(kk))*sin(theta(kk));
      phatj(3)=cos(theta(kk));

      xip=(1-sum(phati.*phatj))/2;
           
      ORF(ll,kk)=3*( 1/3 + xip * ( log(xip) -1/6) ); %Hellings-Downs formula
      if (ll == kk) 
          ORF(ll,kk)=2; %2 along diagonal because we include GW self noise at pulsar 
      end

    end
end

% Use Cholesky transform to take 'square root' of ORF matrix
M=chol(ORF,'lower');

% Calculate frequency dependent pre-factor C(f)
C=H0^2/(16*pi^2)/(2*pi)^2 * f.^(-5) .* Omega * dur;  

% injection residuals in the frequency domain
Res_f=(M * wlocal);
UCRes_f=zeros(Npulsars,Nf);
for ll=1:Npulsars
    Res_f(ll,:) = Res_f(ll,:) .* (C.^(1/2));    %rescale by frequency dependent factor
    Res_f(ll,1)=0;                            %set DC bin to zero to avoid infinities
    Res_f(ll,Nf)=0;                           %set Nyquist bin to zero also (it just needs to be real but lets not add crud to our calculation)
end

% Now fill in bins after Nyquist (for fft data packing) and take inverse FT
% to calculate net residuals
for ll=1:Npulsars
     for kk=1:Nf-2
         Res_f(ll,Nf+kk)=conj(Res_f(ll,Nf-kk));
     end
     
     Res_t(ll,:)=ifft(Res_f(ll,:))/dt;   %ifft includes a factor of 1/N, so divide by dt to effectively multiply by df=1/T
    
end

% add GW residuals to white noise residuals
for ll=1:Npulsars
    Total_Res_t(ll,:)=residualdata(ll,:)+Res_t(ll,:);
end

% write data sets
for ll=1:Npulsars
    pulsarfile=strcat('dataset1-',num2str(ll),'.txt');
    op=Total_Res_t(ll,:);
    save(pulsarfile,'op','-ASCII');
end

% plot total residuals
%figure;plot(ut,Total_Res_t(:,:))



%%%%%%%%%%%%% Generate data with a much smaller background

% Stochastic background parameters in terms of Omega
alpha=0;
Omega_alpha = 5e-12;  %Here I set the value of the background to unity because I want S
Omega=Omega_alpha*f.^alpha;

%create intrinsic white noise residual data for each pulsar    
for ll=1:Npulsars
  residualdata(ll,:)=ResidualRMS*randn(Np,1);    
end

% Create random frequency series from zero mean, unit variance, Gaussian
% distributions
for ll=1:Npulsars
    wlocal(ll,:)=randn(Nf,1)+i*randn(Nf,1);      
end

% Calculate frequency dependent pre-factor C(f)
C=H0^2/(16*pi^2)/(2*pi)^2 * f.^(-5) .* Omega * dur;  

% injection residuals in the frequency domain
Res_f=(M * wlocal);
UCRes_f=zeros(Npulsars,Nf);
for ll=1:Npulsars
    Res_f(ll,:) = Res_f(ll,:) .* (C.^(1/2));    %rescale by frequency dependent factor
    Res_f(ll,1)=0;                            %set DC bin to zero to avoid infinities
    Res_f(ll,Nf)=0;                           %set Nyquist bin to zero also (it just needs to be real but lets not add crud to our calculation)
end

% Now fill in bins after Nyquist (for fft data packing) and take inverse FT
% to calculate net residuals
for ll=1:Npulsars
     for kk=1:Nf-2
         Res_f(ll,Nf+kk)=conj(Res_f(ll,Nf-kk));
     end
     
     Res_t(ll,:)=ifft(Res_f(ll,:))/dt;   %ifft includes a factor of 1/N, so divide by dt to effectively multiply by df=1/T
    
end

% add GW residuals to white noise residuals
for ll=1:Npulsars
    Total_Res_t(ll,:)=residualdata(ll,:)+Res_t(ll,:);
end

% write data sets
for ll=1:Npulsars
    pulsarfile=strcat('dataset2-',num2str(ll),'.txt');
    op=Total_Res_t(ll,:);
    save(pulsarfile,'op','-ASCII');
end

% plot total residuals
%figure;plot(ut,Total_Res_t(:,:))

%%%%%%%%%%%%% Generate uncorrelated data with the same spectrum

% Stochastic background parameters in terms of Omega
alpha=0;
Omega_alpha = 1e-8;  %Here I set the value of the background to unity because I want S
Omega=Omega_alpha*f.^alpha;

%create intrinsic white noise residual data for each pulsar    
for ll=1:Npulsars
  residualdata(ll,:)=ResidualRMS*randn(Np,1);    
end

% Create random frequency series from zero mean, unit variance, Gaussian
% distributions
for ll=1:Npulsars
    wlocal(ll,:)=randn(Nf,1)+i*randn(Nf,1);      
end

% Calculate frequency dependent pre-factor C(f)
C=2*H0^2/(16*pi^2)/(2*pi)^2 * f.^(-5) .* Omega * dur;  

% injection residuals in the frequency domain
Res_f=wlocal;
UCRes_f=zeros(Npulsars,Nf);
for ll=1:Npulsars
    Res_f(ll,:) = Res_f(ll,:) .* (C.^(1/2));    %rescale by frequency dependent factor
    Res_f(ll,1)=0;                            %set DC bin to zero to avoid infinities
    Res_f(ll,Nf)=0;                           %set Nyquist bin to zero also (it just needs to be real but lets not add crud to our calculation)
end

% Now fill in bins after Nyquist (for fft data packing) and take inverse FT
% to calculate net residuals
for ll=1:Npulsars
     for kk=1:Nf-2
         Res_f(ll,Nf+kk)=conj(Res_f(ll,Nf-kk));
     end
     
     Res_t(ll,:)=ifft(Res_f(ll,:))/dt;   %ifft includes a factor of 1/N, so divide by dt to effectively multiply by df=1/T
    
end

% add GW residuals to white noise residuals
for ll=1:Npulsars
    Total_Res_t(ll,:)=residualdata(ll,:)+Res_t(ll,:);
end

% write data sets
for ll=1:Npulsars
    pulsarfile=strcat('dataset3-',num2str(ll),'.txt');
    op=Total_Res_t(ll,:);
    save(pulsarfile,'op','-ASCII');
end

% plot total residuals
%figure;plot(ut,Total_Res_t(:,:))











