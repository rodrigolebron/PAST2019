close all
clear all
clc

%%% Loading Data to Evaluate
Etheta_Ephi_X_Data = csvread('Etheta_Ephi_Xpol.csv',1,0);

THETA = reshape(Etheta_Ephi_X_Data(:,1),181,[]);
PHI = reshape(Etheta_Ephi_X_Data(:,2),181,[]);

Ephi_X = (reshape(Etheta_Ephi_X_Data(:,3),181,[])+1i.*reshape(Etheta_Ephi_X_Data(:,4),181,[]))./1000;
Etheta_X = (reshape(Etheta_Ephi_X_Data(:,5),181,[])+1i.*reshape(Etheta_Ephi_X_Data(:,6),181,[]))./1000;

[u_co_2X,u_cross_2X,u_co_3X,u_cross_3X,u_co_2X_dB,u_cross_2X_dB,u_co_3X_dB,u_cross_3X_dB] = efields2cocrossx(deg2rad(THETA),deg2rad(PHI),Etheta_X,Ephi_X);

vector1X0 = u_co_2X(91,91);
vector2X0 = u_cross_2X(91,91);
vector1X45 = u_co_2X(91+45,91-45);
vector2X45 = u_cross_2X(91+45,91-45);

vector1_magX0 = abs(vector1X0);
vector1_angX0 = angle(vector1X0);

figure
compass(vector1X0/vector1_magX0,'b'); hold on
compass(vector2X0/vector1_magX0,'r')

figure
compass(vector1X45/vector1_magX0,'b'); hold on
compass(vector2X45/vector1_magX0,'r')

figure
compass(vector1X0/vector1_magX0*exp(-1j*vector1_angX0),'b'); hold on
compass(vector2X0/vector1_magX0*exp(-1j*vector1_angX0),'r')

figure
compass(vector1X45/vector1_magX0*exp(-1j*vector1_angX0),'b'); hold on
compass(vector2X45/vector1_magX0*exp(-1j*vector1_angX0),'r')



% figure
% imagesc(real(u_co_2X))
% 
% figure
% imagesc(imag(u_co_2X))
% 
% figure
% imagesc(abs(u_co_2X))
% 
% figure
% imagesc(abs(u_cross_2X))

%%% Plot a Vector and its Components