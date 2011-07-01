%%% This script is to double check the variance of the coherence statistic

clear

load data/SkyPositions.txt

%%%%% Number of pulsars
Npulsars = length(SkyPositions);
Npairs=Npulsars*(Npulsars-1)/2;

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
Omega_alpha = 1e-8;  
Omega=Omega_alpha*f.^alpha;

jj=1;
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

% find expected HD value for each pair (I could be more clever about this I know)
for ll=1:Npulsars 
    
     phati(1)=cos(phi(ll))*sin(theta(ll));
     phati(2)=sin(phi(ll))*sin(theta(ll));
     phati(3)=cos(theta(ll));
    
    for kk=ll+1:Npulsars
      
      phatj(1)=cos(phi(kk))*sin(theta(kk));
      phatj(2)=sin(phi(kk))*sin(theta(kk));
      phatj(3)=cos(theta(kk));
    
      angularseparation(jj)=acos(sum(phati.*phatj));
      xip=(1-sum(phati.*phatj))/2;
      xi(jj)=3*( 1/3 + xip * ( log(xip) -1/6) ); %Hellings-Downs formula

      jj=jj+1;
    end
end

for mm=1:2000
    
    %%%%%% CREATE RESIDUALS %%%%%%%%%%
    %create intrinsic white noise residual data for each pulsar    
    for ll=1:Npulsars
        residualdata(ll,:)=ResidualRMS*randn(Np,1);    
    end
    % Create random frequency series from zero mean, unit variance, Gaussian
    % distributions
    for ll=1:Npulsars
        wlocal(ll,:)=randn(Nf,1)+i*randn(Nf,1);      
    end
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

    %%%%%%% ANALYZE RESIDUALS %%%%%%%%%%
    % find correlation for first dataset
    jj=1;
    for ll=1:Npulsars 
        for kk=ll+1:Npulsars
            correlation(jj)=1/Np*sum(Total_Res_t(ll,:).*Total_Res_t(kk,:));     
            jj=jj+1;
        end
    end

    %%%%%%%% calculate the coherence
    % first calculate the mean values of the correlation and expected correlations
    meanr=mean(correlation);
    meanxi=mean(xi);

    % then calculate the standard deviations
    stdr=sqrt(1/Npairs * sum((correlation-meanr).^2));
    stdxi=sqrt(1/Npairs * sum((xi-meanxi).^2));

    % finally calculate the coherence
    rho(mm)=1/Npairs * sum( (correlation-meanr).*(xi-meanxi) )/(stdr*stdxi);

    % compute the significance of the coherence
    Significance(mm)=rho(mm)*sqrt(Npairs);
    mm

end

meanrho=mean(rho)
meanSig=mean(Significance)
stdrho=std(rho)
expectedstdrho=sqrt(1/Npairs)









