clear

directoryname='data/';
load data/SkyPositions.txt

phi=SkyPositions(:,1);
theta=SkyPositions(:,2);

Npulsars=length(SkyPositions);
Np=500;

jj=1;
% find angular separation, and expected HD value for each pair
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
      xi(jj)=3*( 1/3 + xip * ( log(xip) -1/6) )/2; %Hellings-Downs formula

      jj=jj+1;
    end
end

%%%%%%%%%%%%%%%% FIRST DATA SET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the first pulsar data set
for ll=1:Npulsars
    pulsarfile=strcat(directoryname,'dataset1-',num2str(ll),'.txt');
    residualdata(ll,:)=load(pulsarfile);
end

% plot residuals
figure(1);plot(1:500,residualdata(:,:))

% find correlation for first dataset
jj=1;
for ll=1:Npulsars 
    for kk=ll+1:Npulsars
      correlation(jj)=1/Np*sum(residualdata(ll,:).*residualdata(kk,:));     
      jj=jj+1;
    end
end

% plot the correlation
figure(2); plot(angularseparation,correlation,'.');
Npairs=Npulsars*(Npulsars-1)/2;

%%%%%%%% calculate the coherence
% first calculate the mean values of the correlation and expected correlations
meanr=mean(correlation);
meanxi=mean(xi);

% then calculate the standard deviations
stdr=sqrt(1/Npairs * sum((correlation-meanr).^2));
stdxi=sqrt(1/Npairs * sum((xi-meanxi).^2));

% finally calculate the coherence
rho=1/Npairs * sum( (correlation-meanr).*(xi-meanxi) )/(stdr*stdxi)

% compute the significance of the coherence
Significance=rho*sqrt(Npairs)

%%%%%%%%%%%%%%%% SECOND DATA SET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the first pulsar data set
for ll=1:Npulsars
    pulsarfile=strcat(directoryname,'dataset2-',num2str(ll),'.txt');
    residualdata(ll,:)=load(pulsarfile);
end

% plot residuals
figure(3);plot(1:500,residualdata(:,:))

% find correlation for first dataset
jj=1;
for ll=1:Npulsars 
    for kk=ll+1:Npulsars
      correlation(jj)=1/Np*sum(residualdata(ll,:).*residualdata(kk,:));     
      jj=jj+1;
    end
end

% plot the correlation
figure(4); plot(angularseparation,correlation,'.');
Npairs=Npulsars*(Npulsars-1)/2;

%%%%%%%% calculate the coherence
% first calculate the mean values of the correlation and expected correlations
meanr=mean(correlation);
meanxi=mean(xi);

% then calculate the standard deviations
stdr=sqrt(1/Npairs * sum((correlation-meanr).^2));
stdxi=sqrt(1/Npairs * sum((xi-meanxi).^2));

% finally calculate the coherence
rho=1/Npairs * sum( (correlation-meanr).*(xi-meanxi) )/(stdr*stdxi)

% compute the significance of the coherence
Significance=rho*sqrt(Npairs)


%%%%%%%%%%%%%%%% THIRD DATA SET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the first pulsar data set
for ll=1:Npulsars
    pulsarfile=strcat(directoryname,'dataset3-',num2str(ll),'.txt');
    residualdata(ll,:)=load(pulsarfile);
end

% plot residuals
figure(5);plot(1:500,residualdata(:,:))

% find correlation for first dataset
jj=1;
for ll=1:Npulsars 
    for kk=ll+1:Npulsars
      correlation(jj)=1/Np*sum(residualdata(ll,:).*residualdata(kk,:));     
      jj=jj+1;
    end
end

% plot the correlation
figure(6); plot(angularseparation,correlation,'.');
Npairs=Npulsars*(Npulsars-1)/2;

%%%%%%%% calculate the coherence
% first calculate the mean values of the correlation and expected correlations
meanr=mean(correlation);
meanxi=mean(xi);

% then calculate the standard deviations
stdr=sqrt(1/Npairs * sum((correlation-meanr).^2));
stdxi=sqrt(1/Npairs * sum((xi-meanxi).^2));

% finally calculate the coherence
rho=1/Npairs * sum( (correlation-meanr).*(xi-meanxi) )/(stdr*stdxi)

% compute the significance of the coherence
Significance=rho*sqrt(Npairs)





