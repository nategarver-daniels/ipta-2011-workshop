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