%script that plots antenna pattern response of pulsar-earth system and LIGO

% the idea here is we give the students a script that plots a meshed sphere and ask them to
% edit it to show the pulsar earth system response
clear

phi=linspace(0,2*pi,100);
theta=linspace(0,pi,100);

% set up grid
[phi,theta]=meshgrid(phi,theta);

% radius of sphere
r=5;

x=r.*sin(theta).*cos(phi);
y=r.*sin(theta).*sin(phi);
z=r.*cos(theta);

% plot sphere
figure; mesh(x,y,z)

% they just need to make the radius a function of theta and phi
% antenna pattern response to a +-polarized wave
% Eq. 7 of Anholm et al. (Phys. Rev. D 79, 084030 (2009))
rp=abs(1/2*sin(theta).^2 .*(cos(phi).^2-sin(phi).^2) ./ (1+cos(theta)));

x=rp.*sin(theta).*cos(phi);
y=rp.*sin(theta).*sin(phi);
z=rp.*cos(theta);

figure; mesh(x,y,z)

% antenna pattern including frequency dependence (i.e. the pulsar term)
% (Eqs. 16 and 17 of Anholm et al.)

fL=10; % this is typical for pulsars
%fL=1; % for students to try. They should notice a frequency dependence in
%the picture

rp=abs(    (exp(-2*pi*i*fL*(1+cos(theta)))-1) .*...
    1/2*sin(theta).^2 .*(cos(phi).^2-sin(phi).^2)./ (1+cos(theta)));

x=rp.*sin(theta).*cos(phi);
y=rp.*sin(theta).*sin(phi);
z=rp.*cos(theta);

figure; mesh(x,y,z)

% EXCERCISE: At this point they do a Taylor expansion for theta=pi+delta, where delta
% is sufficiebtly small that the exponential can be expanded. This way they
% see where the frequency dependence in the response at the south pole is coming from 

% antenna pattern for LIGO (see Eq. 3 of Phys. Rev. D 79, 082002 (2009))
% this is kind of a cool paper, because it deals with non-Einstein tehories
% of gravity

% plus polarization with the polarization angle set to 0

psi=0;

rp=abs( 1/2 * (1+cos(theta).^2) .* cos(2*phi) .* cos(2*psi)...
    -cos(theta).*sin(2*theta)*sin(2*psi));

x=rp.*sin(theta).*cos(phi);
y=rp.*sin(theta).*sin(phi);
z=rp.*cos(theta);

figure; mesh(x,y,z)






