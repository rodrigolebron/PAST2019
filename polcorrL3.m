function [ExcitationX,ExcitationY,ETHETA_TOTAL,EPHI_TOTAL] = polcorrL3(THETA,PHI,ETHETA_X,EPHI_X,ETHETA_Y,EPHI_Y,theta0,phi0,Ex_goal,Ey_goal)
%Example
%[ExcitationX,ExcitationY,ETHETA_TOTAL,EPHI_TOTAL] = polcorrL3(THETA,PHI,ETHETA_X,EPHI_X,ETHETA_Y,EPHI_Y,theta0,phi0,Ex_goal,Ey_goal)

% %%% Loading Data of Patterns
% Etheta_Ephi_X_Data = csvread('Etheta_Ephi_Xpol.csv',1,0);
% Etheta_Ephi_Y_Data = csvread('Etheta_Ephi_Ypol.csv',1,0);
% 
% %%% Populating Matrix Data for Processing
% THETA = deg2rad(reshape(Etheta_Ephi_X_Data(:,1),181,[]));
% PHI = deg2rad(reshape(Etheta_Ephi_X_Data(:,2),181,[]));
% 
% EPHI_X = (reshape(Etheta_Ephi_X_Data(:,3),181,[])+1i.*reshape(Etheta_Ephi_X_Data(:,4),181,[]))./1000;
% ETHETA_X = (reshape(Etheta_Ephi_X_Data(:,5),181,[])+1i.*reshape(Etheta_Ephi_X_Data(:,6),181,[]))./1000;
% 
% EPHI_Y = (reshape(Etheta_Ephi_Y_Data(:,3),181,[])+1i.*reshape(Etheta_Ephi_Y_Data(:,4),181,[]))./1000;
% ETHETA_Y = (reshape(Etheta_Ephi_Y_Data(:,5),181,[])+1i.*reshape(Etheta_Ephi_Y_Data(:,6),181,[]))./1000;

[~,~,u_co_3X,u_cross_3X,~,~,~,~] = efields2cocrossx(THETA,PHI,ETHETA_X,EPHI_X);
[~,~,u_co_3Y,u_cross_3Y,~,~,~,~] = efields2cocrossy(THETA,PHI,ETHETA_Y,EPHI_Y);

%%% Function Input Parameters

% theta0 = 45;
% phi0 = 45;
% Ex_goal = 1;
% Ey_goal = 0;

THETA_DEG = rad2deg(THETA);
PHI_DEG = rad2deg(PHI);

%%%%%%%%%%%%%%%%%%%%%%%% Finding the Correction Angle based on the Input Data %%%%%%%%%%%%%%%%%%%%%%%%

[MinIth1,Ith1] = min(abs(THETA_DEG - theta0));
[~,Ith2] = min(abs(THETA_DEG' - theta0));

[MinIph1,Iph1] = min(abs(PHI_DEG - phi0));
[~,Iph2] = min(abs(PHI_DEG' - phi0));

if MinIth1(1,1) == MinIth1(1,2)
    Ith0 = Ith1(1,1);
    th_closest = THETA_DEG(Ith0,Ith0);
%     disp(strcat('The closest theta position is:',{' '},num2str(th_closest),{' '},'degrees'))
else
    Ith0 = Ith2(1,1);
    th_closest = THETA_DEG(Ith0,Ith0);
%     disp(strcat('The closest theta position is:',{' '},num2str(th_closest),{' '},'degrees'))

%     THETA = THETA';
%     PHI = PHI';
%     ETHETA_X = ETHETA_X';
%     EPHI_X = EPHI_X';
%     ETHETA_Y = ETHETA_Y';
%     EPHI_Y = EPHI_Y';    
end

if MinIph1(1,1) == MinIph1(1,2)
    Iph0 = Iph1(1,1);
    ph_closest = PHI_DEG(Iph0,Iph0);
%     disp(strcat('The closest phi position is:',{' '},num2str(ph_closest),{' '},'degrees'))
else
    Iph0 = Iph2(1,1);
    ph_closest = PHI_DEG(Iph0,Iph0);    
%     disp(strcat('The closest phi position is:',{' '},num2str(ph_closest),{' '},'degrees'))

%     THETA = THETA';
%     PHI = PHI';
%     ETHETA_X = ETHETA_X';
%     EPHI_X = EPHI_X';
%     ETHETA_Y = ETHETA_Y';
%     EPHI_Y = EPHI_Y';   
end

%%% Applying Correction

Ex_H = u_co_3X(Ith0,Iph0);
Ey_H = u_cross_3X(Ith0,Iph0);
Ex_V = u_cross_3Y(Ith0,Iph0);
Ey_V = u_co_3Y(Ith0,Iph0);

Tx = [Ex_H,Ex_V; Ey_H,Ey_V];
goal_vec = [Ex_goal; Ey_goal;];

V = Tx\goal_vec;
V_Mag = abs(V) .^2 ;

Mags = V_Mag./sqrt(V_Mag(1).^2+V_Mag(2).^2);
Phases = angle(V);

ExcitationX = Mags(1).^.5.*exp(1j.*Phases(1));
ExcitationY = Mags(2).^.5.*exp(1j.*Phases(2));

ETHETA_TOTAL = ExcitationX.*ETHETA_X + ExcitationY.*ETHETA_Y;
EPHI_TOTAL = ExcitationX.*EPHI_X + ExcitationY.*EPHI_Y;

end