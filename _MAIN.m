close all
clear all
clc

Etheta_Ephi = csvread('Etheta_Ephi_Xpol.csv',1,0);
THETA = deg2rad(reshape(Etheta_Ephi(:,1),181,[]));
PHI = deg2rad(reshape(Etheta_Ephi(:,2),181,[]));
Ephi = (reshape(Etheta_Ephi(:,3),181,[])+1i.*reshape(Etheta_Ephi(:,4),181,[]))./1000;
Etheta = (reshape(Etheta_Ephi(:,5),181,[])+1i.*reshape(Etheta_Ephi(:,6),181,[]))./1000;

% N = 11;
% F(N) = struct('cdata',[],'colormap',[]);

% for m = 1 : N

alpha = 00*pi/180;
beta = 0*pi/180;
gamma = 90*pi/180;
Figures = 1;

[THETA_P,PHI_P,Htheta_P,Hphi_P] = thph2thpphp(alpha,beta,gamma,THETA,PHI,Ephi,Etheta,Figures);

% F(m) = getframe(gcf);
% end

% fig = figure;
% movie(fig,F,1)

A = rad2deg(PHI);
AA = rad2deg(PHI_P);

B = rad2deg(THETA);
BB = rad2deg(THETA_P);