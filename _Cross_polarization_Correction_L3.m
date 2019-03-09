close all
clear all
clc

%%% Loading Data of Patterns
Etheta_Ephi_X_Data = csvread('Etheta_Ephi_Xpol.csv',1,0);
Etheta_Ephi_Y_Data = csvread('Etheta_Ephi_Ypol.csv',1,0);

%%% Populating Matrix Data for Processing
THETA = deg2rad(reshape(Etheta_Ephi_X_Data(:,1),181,[]));
PHI = deg2rad(reshape(Etheta_Ephi_X_Data(:,2),181,[]));

Ephi_X = (reshape(Etheta_Ephi_X_Data(:,3),181,[])+1i.*reshape(Etheta_Ephi_X_Data(:,4),181,[]))./1000;
Etheta_X = (reshape(Etheta_Ephi_X_Data(:,5),181,[])+1i.*reshape(Etheta_Ephi_X_Data(:,6),181,[]))./1000;

Ephi_Y = (reshape(Etheta_Ephi_Y_Data(:,3),181,[])+1i.*reshape(Etheta_Ephi_Y_Data(:,4),181,[]))./1000;
Etheta_Y = (reshape(Etheta_Ephi_Y_Data(:,5),181,[])+1i.*reshape(Etheta_Ephi_Y_Data(:,6),181,[]))./1000;

%%% Defining Coordinate Systems
THETA(THETA==0) = 0.00001;
PHI(PHI==0) = 0.00001;

U = sin(THETA).*cos(PHI);
V = sin(THETA).*sin(PHI);

[AZ,EL] = thph2azel(THETA,PHI);
%%% Calculating Co- and Cross-pol in Ludwig's 3rd Definition

[~,~,u_co_3Y,u_cross_3Y,~,~,u_co_3Y_dB,u_cross_3Y_dB] = efields2cocrossy(THETA,PHI,Etheta_Y,Ephi_Y);
[~,~,u_co_3X,u_cross_3X,~,~,u_co_3X_dB,u_cross_3X_dB] = efields2cocrossx(THETA,PHI,Etheta_X,Ephi_X);

[Plot1] = plot3d(rad2deg(AZ),rad2deg(EL),u_co_3X_dB,'AZ','EL',[-50 0]);
[Plot2] = plot3d(rad2deg(AZ),rad2deg(EL),u_cross_3X_dB,'AZ','EL',[-50 0]);
% [Plot3] = plot3d(U,V,u_co_3X_dB,'u','v',[-50 0]);
% [Plot4] = plot3d(U,V,u_cross_3X_dB,'u','v',[-50 0]);

%%% Cross-polarization Correction Algorithm in Transmit

theta_cal = 45;
phi_cal = 45;

theta_cal = deg2rad(theta_cal);
phi_cal = deg2rad(phi_cal);

[theta_row,theta_col] = find(THETA == theta_cal);
[phi_row,phi_col] = find(PHI == phi_cal);

Ex_H = u_co_3X(theta_row(1,1),phi_col(1,1));
Ey_H = u_cross_3X(theta_row(1,1),phi_col(1,1));
Ex_V = u_cross_3Y(theta_row(1,1),phi_col(1,1));
Ey_V = u_co_3Y(theta_row(1,1),phi_col(1,1));

%u_cross_3Y(theta_row(1,1),phi_col(1,1))
% 1. Solve for the D-plane pattern cuts for realized gain (L3X / L3Y) in HFSS
% 2. make a data table for a far field report AT JUST THE DESIRED CORRECTION ANGLE:  get re(rEL3X), im(rEL3X), re(rEL3Y), im(rEL3Y)
% 3. Go into HFSS to edit the sources, turn on just the H port (the one along X axis)
% 4. Enter the values for Ex_H and Ey_H below
% 5. repeat edit sources, turn on just the V port (along Y axis)
% 6. repeat fill in values for Ex_V and Ey_V below
% 7. save and run the matlab script
% 8. put the output values from the script into the field source in HFSS

% Ex_H = (13.024375)+j*(2.980399);
% Ey_H = (-80.361649/1000)+j*(295.145849/1000);
% Ex_V = (-79.259259/1000)+j*(295.103508/1000);
% Ey_V = (13.025491)+j*(2.975170);

Tx = [Ex_H Ex_V; Ey_H Ey_V; ];
Tx_inv = inv(Tx);

Ex_goal = 1;
Ey_goal = 0;

goal_vec = [Ex_goal; Ey_goal;];

v = Tx_inv*goal_vec;
v_mag = abs(v) .^2 ; % excell spreadsheet and hfss need this value squared

Mags = v_mag./sqrt(v_mag(1).^2+v_mag(2).^2) % these values get put into HFSS source magnitude fields
Phases = angle(v) % these values get put into the HFSS source angle fields

Excitation1 = Mags(1).^.5.*exp(1j.*Phases(1));
Excitation2 = Mags(2).^.5.*exp(1j.*Phases(2));

Etheta_total = Excitation1.*Etheta_X + Excitation2.*Etheta_Y;
Ephi_total = Excitation1.*Ephi_X + Excitation2.*Ephi_Y;

Etheta_total(theta_row(1,1),phi_col(1,1))
Ephi_total(theta_row(1,1),phi_col(1,1))

[~,~,u_co_3X,u_cross_3X,~,~,u_co_3X_dB,u_cross_3X_dB] = efields2cocrossx(THETA,PHI,Etheta_total,Ephi_total);

[Plot5] = plot3d(rad2deg(AZ),rad2deg(EL),u_co_3X_dB,'AZ','EL',[-50 0]);
[Plot6] = plot3d(rad2deg(AZ),rad2deg(EL),u_cross_3X_dB,'AZ','EL',[-50 0]);
% [Plot7] = plot3d(U,V,u_co_3X_dB,'u','v',[-50 0]);
% [Plot8] = plot3d(U,V,u_cross_3X_dB,'u','v',[-50 0]);

[Pattern_Normalized,AF2D] = A2DS(3e9,32,32,50e-3,50e-3,rad2deg(theta_cal),rad2deg(phi_cal),THETA,PHI);

[Plot9] = plot3d(U,V,AF2D,'u','v',[-50 0]);

u_co_3X_normalized = abs(u_co_3X)./max(max(abs(u_co_3X)));
u_cross_3X_normalized = abs(u_cross_3X)./max(max(abs(u_co_3X)));

Array_Pattern_Co = 10.*log10((Pattern_Normalized.*u_co_3X_normalized).^2);
Array_Pattern_Cross = 10.*log10((Pattern_Normalized.*u_cross_3X_normalized).^2);

[Plot10] = plot3d(U,V,Array_Pattern_Co,'u','v',[-50 0]);
[Plot11] = plot3d(U,V,Array_Pattern_Cross,'u','v',[-50 0]);


