close all
clear all
clc

%%% Loading Data of Patterns
Etheta_Ephi_X_Data = csvread('Etheta_Ephi_Xpol.csv',1,0);
Etheta_Ephi_Y_Data = csvread('Etheta_Ephi_Ypol.csv',1,0);

%%% Populating Matrix Data for Processing
THETA = deg2rad(reshape(Etheta_Ephi_X_Data(:,1),181,[]));
PHI = deg2rad(reshape(Etheta_Ephi_X_Data(:,2),181,[]));

[AZ,EL] = thph2azel(THETA,PHI);

EPHI_X = (reshape(Etheta_Ephi_X_Data(:,3),181,[])+1i.*reshape(Etheta_Ephi_X_Data(:,4),181,[]))./1000;
ETHETA_X = (reshape(Etheta_Ephi_X_Data(:,5),181,[])+1i.*reshape(Etheta_Ephi_X_Data(:,6),181,[]))./1000;
EPHI_X_dB = 10.*log10(abs(EPHI_X).^2);
ETHETA_X_dB = 10.*log10(abs(ETHETA_X).^2);

Normalization_X = max(max(EPHI_X_dB));

EPHI_Y = (reshape(Etheta_Ephi_Y_Data(:,3),181,[])+1i.*reshape(Etheta_Ephi_Y_Data(:,4),181,[]))./1000;
ETHETA_Y = (reshape(Etheta_Ephi_Y_Data(:,5),181,[])+1i.*reshape(Etheta_Ephi_Y_Data(:,6),181,[]))./1000;
EPHI_Y_dB = 10.*log10(abs(EPHI_Y).^2);
ETHETA_Y_dB = 10.*log10(abs(ETHETA_Y).^2);

Normalization_Y = max(max(EPHI_Y_dB));

[u_co_2Y,u_cross_2Y,u_co_3Y,u_cross_3Y,u_co_2Y_dB,u_cross_2Y_dB,u_co_3Y_dB,u_cross_3Y_dB] = efields2cocrossy(THETA,PHI,ETHETA_Y,EPHI_Y);
[u_co_2X,u_cross_2X,u_co_3X,u_cross_3X,u_co_2X_dB,u_cross_2X_dB,u_co_3X_dB,u_cross_3X_dB] = efields2cocrossx(THETA,PHI,ETHETA_X,EPHI_X);

% [Plot1] = plot3d(rad2deg(AZ),rad2deg(EL),ETHETA_X_dB-Normalization_X,'Azimuth','Elevation',[-40 0]);
% [Plot2] = plot3d(rad2deg(AZ),rad2deg(EL),EPHI_X_dB-Normalization_X,'Azimuth','Elevation',[-40 0]);
% [Plot3] = plot3d(rad2deg(AZ),rad2deg(EL),u_co_2X_dB,'Azimuth','Elevation',[-40 0]);
% [Plot4] = plot3d(rad2deg(AZ),rad2deg(EL),u_cross_2X_dB,'Azimuth','Elevation',[-40 0]);
% [Plot5] = plot3d(rad2deg(AZ),rad2deg(EL),u_co_3X_dB,'Azimuth','Elevation',[-40 0]);
% [Plot6] = plot3d(rad2deg(AZ),rad2deg(EL),u_cross_3X_dB,'Azimuth','Elevation',[-40 0]);

% [Plot7] = plot3d(rad2deg(AZ),rad2deg(EL),ETHETA_Y_dB-Normalization_Y,'Azimuth','Elevation',[-40 0]);
% [Plot8] = plot3d(rad2deg(AZ),rad2deg(EL),EPHI_Y_dB-Normalization_Y,'Azimuth','Elevation',[-40 0]);
% [Plot9] = plot3d(rad2deg(AZ),rad2deg(EL),u_co_2Y_dB,'Azimuth','Elevation',[-40 0]);
% [Plot10] = plot3d(rad2deg(AZ),rad2deg(EL),u_cross_2Y_dB,'Azimuth','Elevation',[-40 0]);
% [Plot11] = plot3d(rad2deg(AZ),rad2deg(EL),u_co_3Y_dB,'Azimuth','Elevation',[-40 0]);
% [Plot12] = plot3d(rad2deg(AZ),rad2deg(EL),u_cross_3Y_dB,'Azimuth','Elevation',[-40 0]);


Condition_Number_Matrix_ThPh = zeros(181,181);
Condition_Number_Matrix_2L = zeros(181,181);
Condition_Number_Matrix_3L = zeros(181,181);


for m = 1 : 181
    for n = 1 : 181

        Jones_Matrix_ThPh = [EPHI_X(m,n),EPHI_Y(m,n);ETHETA_X(m,n),ETHETA_Y(m,n)];
        Condition_Number_Matrix_ThPh(m,n) = cond(Jones_Matrix_ThPh);  

%         Jones_Matrix_2L = [u_co_2X(m,n),u_cross_2Y(m,n);u_cross_2X(m,n),u_co_2Y(m,n)];
%         Condition_Number_Matrix_2L(m,n) = cond(Jones_Matrix_2L);  
        
        Jones_Matrix_3L = [u_co_3X(m,n),u_cross_3Y(m,n);u_cross_3X(m,n),u_co_3Y(m,n)];
        Condition_Number_Matrix_3L(m,n) = cond(Jones_Matrix_3L);       
        
    end
end

IXR_ThPh = ((Condition_Number_Matrix_ThPh+1)./(Condition_Number_Matrix_ThPh-1)).^2;
% IXR_2L = ((Condition_Number_Matrix_2L+1)./(Condition_Number_Matrix_2L-1)).^2;
IXR_3L = ((Condition_Number_Matrix_3L+1)./(Condition_Number_Matrix_3L-1)).^2;

IXR_ThPh_dB = 10.*log10(IXR_ThPh);
IXR_Norm_dB_ThPh = max(max(IXR_ThPh_dB));

% IXR_2L_dB = 10.*log10(IXR_2L);
% IXR_Norm_dB_2L = max(max(IXR_2L_dB));

IXR_3L_dB = 10.*log10(IXR_3L);
IXR_Norm_dB_3L = max(max(IXR_3L_dB));

[Plot13] = plot3d(rad2deg(AZ),rad2deg(EL),IXR_ThPh_dB,'Azimuth','Elevation',[0 60]);
% [Plot14] = plot3d(rad2deg(AZ),rad2deg(EL),IXR_2L_dB,'Azimuth','Elevation',[0 60]);
[Plot15] = plot3d(rad2deg(AZ),rad2deg(EL),IXR_3L_dB,'Azimuth','Elevation',[0 60]);