function [] = pol_corr(theta0,phi0,Ex_goal,Ey_goal,THETA,PHI,ETHETA_X,EPHI_X,ETHETA_Y,EPHI_Y,Figures)

THETA_DEG = rad2deg(THETA);
PHI_DEG = rad2deg(PHI);

[AZ,EL] = thph2azel(THETA,PHI);

%%%%%%%%%%%%%%%%%%%%%%%% Finding the Correction Angle based on the Input Data %%%%%%%%%%%%%%%%%%%%%%%%

[MinIth1,Ith1] = min(abs(THETA_DEG - theta0));
[~,Ith2] = min(abs(THETA_DEG' - theta0));

[MinIph1,Iph1] = min(abs(PHI_DEG - phi0));
[~,Iph2] = min(abs(PHI_DEG' - phi0));

if MinIth1(1,1) == MinIth1(1,2)
    Ith0 = Ith1(1,1);
    th_closest = THETA_DEG(Ith0,Ith0);
    disp(strcat('The closest theta position is:',{' '},num2str(th_closest),{' '},'degrees'))
else
    Ith0 = Ith2(1,1);
    th_closest = THETA_DEG(Ith0,Ith0);
    disp(strcat('The closest theta position is:',{' '},num2str(th_closest),{' '},'degrees'))
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
    disp(strcat('The closest phi position is:',{' '},num2str(ph_closest),{' '},'degrees'))
else
    Iph0 = Iph2(1,1);
    ph_closest = PHI_DEG(Iph0,Iph0);    
    disp(strcat('The closest phi position is:',{' '},num2str(ph_closest),{' '},'degrees'))
%     THETA = THETA';
%     PHI = PHI';
%     ETHETA_X = ETHETA_X';
%     EPHI_X = EPHI_X';
%     ETHETA_Y = ETHETA_Y';
%     EPHI_Y = EPHI_Y';   
end

%%%%%%%%%%%%%%%%%%%%%%%% Calculating Array Factor %%%%%%%%%%%%%%%%%%%%%%%%

[Array_Pattern,~] = A2DS(3e9,32,32,50e-3,50e-3,th_closest,ph_closest,THETA_DEG,PHI_DEG);

%%%%%%%%%%%%%%%%%%%%%%%% Calculating Polarization Basis %%%%%%%%%%%%%%%%%%%%%%%%

[u_co_2Y_uncal,u_cross_2Y_uncal,u_co_3Y_uncal,u_cross_3Y_uncal,u_co_2Y_dB_uncal,u_cross_2Y_dB_uncal,u_co_3Y_dB_uncal,u_cross_3Y_dB_uncal] = efields2cocrossy(THETA,PHI,ETHETA_Y,EPHI_Y);
[u_co_2X_uncal,u_cross_2X_uncal,u_co_3X_uncal,u_cross_3X_uncal,u_co_2X_dB_uncal,u_cross_2X_dB_uncal,u_co_3X_dB_uncal,u_cross_3X_dB_uncal] = efields2cocrossx(THETA,PHI,ETHETA_X,EPHI_X);

%%%%%%%%%%%%%%%%%%%%%%%% Algorithm Starts %%%%%%%%%%%%%%%%%%%%%%%%

%% Correction Angle Position Introduced to Matrix Form in Ludwig's 2nd Definition
Ex_H_2 = u_co_2X_uncal(Ith0,Iph0);
Ey_H_2 = u_cross_2X_uncal(Ith0,Iph0);
Ex_V_2 = u_cross_2Y_uncal(Ith0,Iph0);
Ey_V_2 = u_co_2Y_uncal(Ith0,Iph0);

Tx_2 = [Ex_H_2,Ex_V_2; Ey_H_2,Ey_V_2];
Tx_inv_2 = inv(Tx_2);
goal_vec_2 = [Ex_goal; Ey_goal;];

v_2 = Tx_inv_2*goal_vec_2;
v_mag_2 = abs(v_2) .^2 ;

Mags_2 = v_mag_2./sqrt(v_mag_2(1).^2+v_mag_2(2).^2);
Phases_2 = angle(v_2);

Excitation1_2 = Mags_2(1).^.5.*exp(1j.*Phases_2(1));
Excitation2_2 = Mags_2(2).^.5.*exp(1j.*Phases_2(2));

Etheta_total_2 = Excitation1_2.*ETHETA_X + Excitation2_2.*ETHETA_Y;
Ephi_total_2 = Excitation1_2.*EPHI_X + Excitation2_2.*EPHI_Y;

[u_co_2Y_cal,u_cross_2Y_cal,~,~,u_co_2Y_cal_dB,u_cross_2Y_cal_dB,~,~] = efields2cocrossy(THETA,PHI,Etheta_total_2,Ephi_total_2);
[u_co_2X_cal,u_cross_2X_cal,~,~,u_co_2X_cal_dB,u_cross_2X_cal_dB,~,~] = efields2cocrossx(THETA,PHI,Etheta_total_2,Ephi_total_2);

%% Applying Array Factor to Un-calibrated Patterns in Ludwig's 2nd Definition
u_co_2X_uncal_array_pattern = 10.*log10(abs(Array_Pattern.*u_co_2X_uncal).^2);
u_cross_2X_uncal_array_pattern = 10.*log10(abs(Array_Pattern.*u_cross_2X_uncal).^2);
u_co_2Y_uncal_array_pattern = 10.*log10(abs(Array_Pattern.*u_co_2Y_uncal).^2);
u_cross_2Y_uncal_array_pattern = 10.*log10(abs(Array_Pattern.*u_cross_2Y_uncal).^2);

u_co_2X_uncal_array_pattern(u_co_2X_uncal_array_pattern==Inf) = -50;
u_cross_2X_uncal_array_pattern(u_cross_2X_uncal_array_pattern==Inf) = -50;
u_co_2Y_uncal_array_pattern(u_co_2Y_uncal_array_pattern==Inf) = -50;
u_cross_2Y_uncal_array_pattern(u_cross_2Y_uncal_array_pattern==Inf) = -50;

Normalization_uncal_1 = max(max(u_co_2X_uncal_array_pattern));
Normalization_uncal_2 = max(max(u_cross_2X_uncal_array_pattern));
Normalization_uncal_3 = max(max(u_co_2Y_uncal_array_pattern));
Normalization_uncal_4 = max(max(u_cross_2Y_uncal_array_pattern));

Normalization_uncal_L2 = max([Normalization_uncal_1,Normalization_uncal_2,Normalization_uncal_3,Normalization_uncal_4]);

%% Applying Array Factor to Calibrated Patterns in Ludwig's 2nd Definition
u_co_2X_cal_array_pattern = 10.*log10(abs(Array_Pattern.*u_co_2X_cal).^2);
u_cross_2X_cal_array_pattern = 10.*log10(abs(Array_Pattern.*u_cross_2X_cal).^2);
u_co_2Y_cal_array_pattern = 10.*log10(abs(Array_Pattern.*u_co_2Y_cal).^2);
u_cross_2Y_cal_array_pattern = 10.*log10(abs(Array_Pattern.*u_cross_2Y_cal).^2);

u_co_2X_cal_array_pattern(u_co_2X_cal_array_pattern==Inf) = -50;
u_cross_2X_cal_array_pattern(u_cross_2X_cal_array_pattern==Inf) = -50;
u_co_2Y_cal_array_pattern(u_co_2Y_cal_array_pattern==Inf) = -50;
u_cross_2Y_cal_array_pattern(u_cross_2Y_cal_array_pattern==Inf) = -50;

Normalization_cal_1 = max(max(u_co_2X_uncal_array_pattern));
Normalization_cal_2 = max(max(u_cross_2X_uncal_array_pattern));
Normalization_cal_3 = max(max(u_co_2Y_uncal_array_pattern));
Normalization_cal_4 = max(max(u_cross_2Y_uncal_array_pattern));

Normalization_cal_L2 = max([Normalization_cal_1,Normalization_cal_2,Normalization_cal_3,Normalization_cal_4]);

%% Correction Angle Position Introduced to Matrix Form in Ludwig's 3rd Definition
Ex_H_3 = u_co_3X_uncal(Ith0,Iph0);
Ey_H_3 = u_cross_3X_uncal(Ith0,Iph0);
Ex_V_3 = u_cross_3Y_uncal(Ith0,Iph0);
Ey_V_3 = u_co_3Y_uncal(Ith0,Iph0);

Tx_3 = [Ex_H_3,Ex_V_3; Ey_H_3,Ey_V_3];
Tx_inv_3 = inv(Tx_3);
goal_vec_3 = [Ex_goal; Ey_goal;];

v_3 = Tx_inv_3*goal_vec_3;
v_mag_3 = abs(v_3) .^2 ;

Mags_3 = v_mag_3./sqrt(v_mag_3(1).^2+v_mag_3(2).^2);
Phases_3 = angle(v_3);

Excitation1_3 = Mags_3(1).^.5.*exp(1j.*Phases_3(1));
Excitation2_3 = Mags_3(2).^.5.*exp(1j.*Phases_3(2));

Etheta_total_3 = Excitation1_3.*ETHETA_X + Excitation2_3.*ETHETA_Y;
Ephi_total_3 = Excitation1_3.*EPHI_X + Excitation2_3.*EPHI_Y;

[~,~,u_co_3Y_cal,u_cross_3Y_cal,~,~,u_co_3Y_cal_dB,u_cross_3Y_cal_dB] = efields2cocrossy(THETA,PHI,Etheta_total_3,Ephi_total_3);
[~,~,u_co_3X_cal,u_cross_3X_cal,~,~,u_co_3X_cal_dB,u_cross_3X_cal_dB] = efields2cocrossx(THETA,PHI,Etheta_total_3,Ephi_total_3);

%% Applying Array Factor to Un-calibrated Patterns in Ludwig's 3rd Definition
u_co_3X_uncal_array_pattern = 10.*log10(abs(Array_Pattern.*u_co_3X_uncal).^2);
u_cross_3X_uncal_array_pattern = 10.*log10(abs(Array_Pattern.*u_cross_3X_uncal).^2);
u_co_3Y_uncal_array_pattern = 10.*log10(abs(Array_Pattern.*u_co_3Y_uncal).^2);
u_cross_3Y_uncal_array_pattern = 10.*log10(abs(Array_Pattern.*u_cross_3Y_uncal).^2);

Normalization_uncal_1 = max(max(u_co_3X_uncal_array_pattern));
Normalization_uncal_2 = max(max(u_cross_3X_uncal_array_pattern));
Normalization_uncal_3 = max(max(u_co_3Y_uncal_array_pattern));
Normalization_uncal_4 = max(max(u_cross_3Y_uncal_array_pattern));

Normalization_uncal = max([Normalization_uncal_1,Normalization_uncal_2,Normalization_uncal_3,Normalization_uncal_4]);

%% Applying Array Factor to Calibrated Patterns in Ludwig's 3rd Definition
u_co_3X_cal_array_pattern = 10.*log10(abs(Array_Pattern.*u_co_3X_cal).^2);
u_cross_3X_cal_array_pattern = 10.*log10(abs(Array_Pattern.*u_cross_3X_cal).^2);
u_co_3Y_cal_array_pattern = 10.*log10(abs(Array_Pattern.*u_co_3Y_cal).^2);
u_cross_3Y_cal_array_pattern = 10.*log10(abs(Array_Pattern.*u_cross_3Y_cal).^2);

Normalization_cal_1 = max(max(u_co_3X_cal_array_pattern));
Normalization_cal_2 = max(max(u_cross_3X_cal_array_pattern));
Normalization_cal_3 = max(max(u_co_3Y_cal_array_pattern));
Normalization_cal_4 = max(max(u_cross_3Y_cal_array_pattern));

Normalization_cal = max([Normalization_cal_1,Normalization_cal_2,Normalization_cal_3,Normalization_cal_4]);

%% Plotting

if Figures == 1
    
    disp({'Plots Enabled'})
      
    figure
    %%
    %%%subplot1
    subplot(4,4,1)
    surf(rad2deg(AZ),rad2deg(EL),u_co_3X_dB_uncal);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Un-calibrated L3X Co-pol')
    
    %%% subplot2
    subplot(4,4,2)
    surf(rad2deg(AZ),rad2deg(EL),u_cross_3X_dB_uncal);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Un-calibrated L3X Cross-pol')

    %%% subplot3
    subplot(4,4,3)
    surf(rad2deg(AZ),rad2deg(EL),u_co_3Y_dB_uncal);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Un-calibrated L3Y Co-pol')

    %%% subplot4
    subplot(4,4,4)
    surf(rad2deg(AZ),rad2deg(EL),u_cross_3Y_dB_uncal);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Un-calibrated L3Y Cross-pol')
    
    %%% subplot5
    subplot(4,4,5)
    surf(rad2deg(AZ),rad2deg(EL),u_co_3X_uncal_array_pattern-Normalization_uncal);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Un-calibrated L3X Array Co-pol')

    %%% subplot6
    subplot(4,4,6)
    surf(rad2deg(AZ),rad2deg(EL),u_cross_3X_uncal_array_pattern-Normalization_uncal);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Un-calibrated L3X Array Cross-pol') 

    %%% subplot7
    subplot(4,4,7)
    surf(rad2deg(AZ),rad2deg(EL),u_co_3Y_uncal_array_pattern-Normalization_uncal);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Un-calibrated L3Y Array Co-pol')

    %%% subplot8
    subplot(4,4,8)
    surf(rad2deg(AZ),rad2deg(EL),u_cross_3Y_uncal_array_pattern-Normalization_uncal);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Un-calibrated L3Y Array Cross-pol') 
    
    %%% subplot9
    subplot(4,4,9)
    surf(rad2deg(AZ),rad2deg(EL),u_co_3X_cal_dB);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Calibrated L3X Co-pol')
    
    %%% subplot10
    subplot(4,4,10)
    surf(rad2deg(AZ),rad2deg(EL),u_cross_3X_cal_dB);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Calibrated L3X Cross-pol')
    
    %%% subplot11
    subplot(4,4,11)
    surf(rad2deg(AZ),rad2deg(EL),u_co_3Y_cal_dB);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Calibrated L3Y Co-pol')
    
    %%% subplot12
    subplot(4,4,12)
    surf(rad2deg(AZ),rad2deg(EL),u_cross_3Y_cal_dB);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Calibrated L3Y Cross-pol')
    
    %%% subplot13
    subplot(4,4,13)
    surf(rad2deg(AZ),rad2deg(EL),u_co_3X_cal_array_pattern-Normalization_cal);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Calibrated L3X Array Co-pol')

    %%% subplot14
    subplot(4,4,14)
    surf(rad2deg(AZ),rad2deg(EL),u_cross_3X_cal_array_pattern-Normalization_cal);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Calibrated L3X Array Cross-pol') 

    %%% subplot15
    subplot(4,4,15)
    surf(rad2deg(AZ),rad2deg(EL),u_co_3Y_cal_array_pattern-Normalization_cal);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Calibrated L3Y Array Co-pol')

    %%% subplot16
    subplot(4,4,16)
    surf(rad2deg(AZ),rad2deg(EL),u_cross_3Y_cal_array_pattern-Normalization_cal);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Calibrated L3Y Array Cross-pol')     
    %%
    
    figure
    %%
    %%%subplot1
    subplot(4,4,1)
    surf(rad2deg(AZ),rad2deg(EL),u_co_2X_dB_uncal);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Un-calibrated L2X Co-pol')
    
    %%% subplot2
    subplot(4,4,2)
    surf(rad2deg(AZ),rad2deg(EL),u_cross_2X_dB_uncal);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Un-calibrated L2X Cross-pol')

    %%% subplot3
    subplot(4,4,3)
    surf(rad2deg(AZ),rad2deg(EL),u_co_2Y_dB_uncal);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Un-calibrated L2Y Co-pol')

    %%% subplot4
    subplot(4,4,4)
    surf(rad2deg(AZ),rad2deg(EL),u_cross_2Y_dB_uncal);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Un-calibrated L2Y Cross-pol')
    
    %%% subplot5
    subplot(4,4,5)
    surf(rad2deg(AZ),rad2deg(EL),u_co_2X_uncal_array_pattern-Normalization_uncal_L2);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Un-calibrated L2X Array Co-pol')

    %%% subplot6
    subplot(4,4,6)
    surf(rad2deg(AZ),rad2deg(EL),u_cross_2X_uncal_array_pattern-Normalization_uncal_L2);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Un-calibrated L2X Array Cross-pol') 

    %%% subplot7
    subplot(4,4,7)
    surf(rad2deg(AZ),rad2deg(EL),u_co_2Y_uncal_array_pattern-Normalization_uncal_L2);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Un-calibrated L2Y Array Co-pol')

    %%% subplot8
    subplot(4,4,8)
    surf(rad2deg(AZ),rad2deg(EL),u_cross_2Y_uncal_array_pattern-Normalization_uncal_L2);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Un-calibrated L2Y Array Cross-pol') 
    
    %%% subplot9
    subplot(4,4,9)
    surf(rad2deg(AZ),rad2deg(EL),u_co_2X_cal_dB);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Calibrated L2X Co-pol')
    
    %%% subplot10
    subplot(4,4,10)
    surf(rad2deg(AZ),rad2deg(EL),u_cross_2X_cal_dB);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Calibrated L2X Cross-pol')
    
    %%% subplot11
    subplot(4,4,11)
    surf(rad2deg(AZ),rad2deg(EL),u_co_2Y_cal_dB);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Calibrated L2Y Co-pol')
    
    %%% subplot12
    subplot(4,4,12)
    surf(rad2deg(AZ),rad2deg(EL),u_cross_2Y_cal_dB);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Calibrated L2Y Cross-pol')
    
    %%% subplot13
    subplot(4,4,13)
    surf(rad2deg(AZ),rad2deg(EL),u_co_2X_cal_array_pattern-Normalization_cal_L2);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Calibrated L2X Array Co-pol')

    %%% subplot14
    subplot(4,4,14)
    surf(rad2deg(AZ),rad2deg(EL),u_cross_2X_cal_array_pattern-Normalization_cal_L2);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Calibrated L2X Array Cross-pol') 

    %%% subplot15
    subplot(4,4,15)
    surf(rad2deg(AZ),rad2deg(EL),u_co_2Y_cal_array_pattern-Normalization_cal_L2);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Calibrated L2Y Array Co-pol')

    %%% subplot16
    subplot(4,4,16)
    surf(rad2deg(AZ),rad2deg(EL),u_cross_2Y_cal_array_pattern-Normalization_cal_L2);
    view([0 90])
    shading interp
    colormap jet
    ylabel('Elevation (deg)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([min(min(rad2deg(AZ))) max(max(rad2deg(AZ)))])
    ylim([min(min(rad2deg(EL))) max(max(rad2deg(EL)))])
    caxis([-40 0])
    set(gca,'fontsize',12)
    colorbar;
    zlim([-40 0])
    axis square
    title('Calibrated L2Y Array Cross-pol')     
    %%
    
else
    
    disp({'Plots Disabled'})
    
end

end
