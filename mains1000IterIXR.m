% Randomly rotated array fast computation code 
% imports patterns from X and Y polarized dipoles
% rotates the patterns for rot_angles=[-90:90];
% create phase shifts for all X,Y positions
% beam forming for target position angle, and co/cross ratio
% plot excitations
% plots result
% steering without using for loops


close all
clear all
clc
%% Program parameters
%constants
THETA_COMP = 2;
PHI_COMP = 1;
X_POL = 1;
Y_POL = 2;

%Configuration Parameters
%Array dimension
NBR_COLS = 32;
NBR_ROWS = 32;
%frequency
FREQ = 3e9;
%Element rotation
minAngle=-90;
maxAngle=90;
%Target beam - in AZ/EL
polarization = 'L3'; %L2,L3
targetAZ = -0;
targetEL = 0;
targetEco2EcrssIndB = 50;

lambda = 3e8/FREQ;
K = 2*pi/lambda;

%% Load Data
Etheta_Ephi_X = csvread('Etheta_Ephi_Xpol.csv',1,0);
Etheta_Ephi_Y = csvread('Etheta_Ephi_Ypol.csv',1,0);
%angles
THETA = deg2rad(reshape(Etheta_Ephi_X(:,1),181,[]));
%THETA = THETA(:,1:180);
PHI = deg2rad(reshape(Etheta_Ephi_X(:,2),181,[]));
%PHI = PHI(:,1:180);
%angles Az, El, U,V
[AZ,EL] = thph2azel(THETA,PHI);
[U,V] = thph2uv(THETA,PHI);
%Fields
EPHI_X = (reshape(Etheta_Ephi_X(:,3),181,[])+1i.*reshape(Etheta_Ephi_X(:,4),181,[]))./1000;
%EPHI_X = EPHI_X(:,1:180);
ETHETA_X = (reshape(Etheta_Ephi_X(:,5),181,[])+1i.*reshape(Etheta_Ephi_X(:,6),181,[]))./1000;
%ETHETA_X = ETHETA_X(:,1:180);

EPHI_Y = (reshape(Etheta_Ephi_Y(:,3),181,[])+1i.*reshape(Etheta_Ephi_Y(:,4),181,[]))./1000;
%EPHI_Y = EPHI_Y(:,1:180);
ETHETA_Y = (reshape(Etheta_Ephi_Y(:,5),181,[])+1i.*reshape(Etheta_Ephi_Y(:,6),181,[]))./1000;
%ETHETA_Y = ETHETA_Y(:,1:180);

Figures = 0;

%IXRdB_Old = IXRthph(THETA,PHI,ETHETA_X,EPHI_X,ETHETA_Y,EPHI_Y,Figures);
%IXRdB_New = IXRthph(THETA_PX,PHI_PX,ETHETA_PX,EPHI_PX,ETHETA_PY,EPHI_PY,Figures);
%% Create array
%coordinates
X_coord_elements = repmat((-lambda/2*floor(NBR_COLS/2)):lambda/2:(lambda/2*floor(NBR_COLS/2)),NBR_ROWS,1);
Y_coord_elements = repmat([(-lambda/2*floor(NBR_ROWS/2)):lambda/2:(lambda/2*floor(NBR_ROWS/2))]',1,NBR_COLS);
nbr_loop=1
for iter =1:nbr_loop
    tic
    %rng(0,'twister');
    Angle_matrix(:,:,iter) = randi([minAngle maxAngle],NBR_ROWS,NBR_COLS);
    %Angle = 0;
    
    Excitations = zeros(NBR_COLS*NBR_ROWS,4);
    
    for i=1:NBR_COLS
        for j=1:NBR_ROWS
            %rotate element patterns
            [ETHETA_X_rot,EPHI_X_rot,ETHETA_Y_rot,EPHI_Y_rot] = rotZthph(Angle_matrix(j,i,iter),ETHETA_X,EPHI_X,ETHETA_Y,EPHI_Y);
            
            %add phase shift in postion
            phsShift = exp(-1i*K*(U*X_coord_elements(j,i)+V*Y_coord_elements(j,i)));
            ETHETA_X_on_ARRAY = ETHETA_X_rot.*phsShift;
            EPHI_X_on_ARRAY = EPHI_X_rot.*phsShift;
            ETHETA_Y_on_ARRAY = ETHETA_Y_rot.*phsShift;
            EPHI_Y_on_ARRAY = EPHI_Y_rot.*phsShift;
            
            
            %store patterns
            Patterns(:,:,i+(j-1)*NBR_COLS,THETA_COMP,X_POL) = ETHETA_X_on_ARRAY; % 1-Th,2-Ph,3-Elmnt,4-Th/PH,5-X/Y
            Patterns(:,:,i+(j-1)*NBR_COLS,PHI_COMP,X_POL) = EPHI_X_on_ARRAY;
            Patterns(:,:,i+(j-1)*NBR_COLS,THETA_COMP,Y_POL) = ETHETA_Y_on_ARRAY;
            Patterns(:,:,i+(j-1)*NBR_COLS,PHI_COMP,Y_POL) = EPHI_Y_on_ARRAY;
            
        end
    end
    
    %  IXR
    [N_theta,N_phi,N,~,~] = size(Patterns);
    Condition_Number_Matrix_ThPh = zeros(N_theta,N_phi);
    Jones_Matrix_Array = zeros(2,2);
    %this loops takes ~18 seconds
    for m = 1 : N_theta
        for n = 1 : N_phi
            %get Jones matrix  Jones=[EPH_X EPH_Y;ETH_X ETH_Y], third dim is elements
            jonesMatrix = squeeze(permute(Patterns(m,n,:,:,:),[1 2 4 5 3])); %1-PH/TH 2-X/Y 3-elements
            jonesMatrixArray = reshape(permute(jonesMatrix,[1 3 2]),N*2,2);
            
            Condition_Number_Matrix_ThPh(m,n) = cond(jonesMatrixArray);
        end
    end
    
    IXR_ThPh(:,:,iter) = ((Condition_Number_Matrix_ThPh+1)./(Condition_Number_Matrix_ThPh-1)).^2;
    IXR_ThPh_dB = 10.*log10(IXR_ThPh(:,:,iter));
    
    %plot3d(rad2deg(THETA),rad2deg(PHI),IXR_ThPh_dB(:,:,iter),'THETA','PHI',[0 60])
    %plot3d(rad2deg(THETA),rad2deg(PHI),IXR_ThPh_dB(:,:,iter),'THETA','PHI',[0 40])
    
    %plot3d(rad2deg(AZ),rad2deg(EL),IXR_ThPh_dB(:,:,iter),'THETA','PHI',[0 60])
    f=plot3d(rad2deg(AZ),rad2deg(EL),IXR_ThPh_dB,'THETA','PHI',[0 40],0);
    filename=strcat('plots/IXRiter',num2str(iter),'.png');
    print(filename,'-dpng')
    close(f)
    toc
    iter/nbr_loop
end
% 
%%% rotate test
save('Data1000iter.mat')


%% STEERING
% get target angles indexes
[targetTHdeg,targetPHdeg] = AzEl2ThPh(targetAZ,targetEL);
[~, TH_index] = min(abs(THETA(:,1) - targetTHdeg*pi/180));
[~, PH_index] = min(abs(PHI(1,:) - targetPHdeg*pi/180));
PHI(TH_index,PH_index)
THETA(TH_index,PH_index)

% %get Jones matrix  Jones=[EPH_X EPH_Y;ETH_X ETH_Y], third dim is elements
%jonesMatrix = squeeze(permute(Patterns(TH_index,PH_index,:,:,:),[1 2 4 5 3])); %1-PH/TH 2-X/Y 3-elements
%jonesMatrixArray = reshape(permute(jonesMatrix,[1 3 2]),N*2,2);
%get inverse Jones matrix  Jones=[ETH_Y -EPH_Y;-ETH_X EPH_X]/determinant, third dim is elements
[N_theta,N_phi,N,~,~] = size(Patterns)
X_indexes = 1:2:N*2;
Y_indexes = 2:2:N*2;
%  lets build the invMatrix
determ = squeeze(Patterns(TH_index,PH_index,:,THETA_COMP,Y_POL))... %ETH_Y
    .*squeeze(Patterns(TH_index,PH_index,:,PHI_COMP,X_POL))...   %EPH_X
     -squeeze(Patterns(TH_index,PH_index,:,PHI_COMP,Y_POL))... %EPH_Y
    .*squeeze(Patterns(TH_index,PH_index,:,THETA_COMP,X_POL)); %ETH_X
invJonesArray(X_indexes,1) = squeeze(Patterns(TH_index,PH_index,:,THETA_COMP,Y_POL))./determ; %ETH_Y
invJonesArray(X_indexes,2) = -squeeze(Patterns(TH_index,PH_index,:,PHI_COMP,Y_POL))./determ; %-EPH_Y
invJonesArray(Y_indexes,1) = -squeeze(Patterns(TH_index,PH_index,:,THETA_COMP,X_POL))./determ; %-ETH_X
invJonesArray(Y_indexes,2) = squeeze(Patterns(TH_index,PH_index,:,PHI_COMP,X_POL))./determ; %EPH_X
%build target E field [EtargetX]
% the final field will be a sum of multiple fields, and assume make Eph=0
targetEco2Ecrsslin  = db2mag(targetEco2EcrssIndB);

if polarization == 'L3'
    target_EL3co = targetEco2Ecrsslin;
    target_EL3x = 1;
    [target_Eth, target_Eph] = L32ThPh(THETA(TH_index,PH_index)*180/pi,...
        PHI(TH_index,PH_index)*180/pi,...
        target_EL3co,target_EL3x);
elseif polarization == 'L2'
    target_EL2Az = targetEco2Ecrsslin;
    target_EL2El = 1;
    [target_Eth, target_Eph] = L2I2ThPh(THETA(TH_index,PH_index)*180/pi,...
        PHI(TH_index,PH_index)*180/pi,...
        target_EL2Az,target_EL2El);
end
EtargetArray = [target_Eph; target_Eth];
%calculate excitations
ExcitationsXYArray =  invJonesArray*EtargetArray;
%rearrange ExcitationsXYArray = [ExcX ExcY]
ExcitationsCmplx =  [ExcitationsXYArray(X_indexes) ExcitationsXYArray(Y_indexes)];

%to multiply excitations without FOR LOOP (dim X/Y 2 to pos 5, dim elemnt 1 to pos 3)
exc = permute(ExcitationsCmplx,[5 4 1 3 2]);
%% Sum patterns
SteeredPatterns = Patterns.*exc;
%Sum
ETHETA_X_TOTAL = squeeze(sum(SteeredPatterns(:,:,:,THETA_COMP,X_POL),3));
EPHI_X_TOTAL = squeeze(sum(SteeredPatterns(:,:,:,PHI_COMP,X_POL),3));
ETHETA_Y_TOTAL = squeeze(sum(SteeredPatterns(:,:,:,THETA_COMP,Y_POL),3));
EPHI_Y_TOTAL = squeeze(sum(SteeredPatterns(:,:,:,PHI_COMP,Y_POL),3));

%sum
ETHETA_TOTAL = squeeze(sum(sum(SteeredPatterns(:,:,:,THETA_COMP,:),3),5));%ETHETA_X_TOTAL+ETHETA_Y_TOTAL;
EPHI_TOTAL = squeeze(sum(sum(SteeredPatterns(:,:,:,PHI_COMP,:),3),5));%EPHI_X_TOTAL+EPHI_Y_TOTAL;

% Eth_sum_X = sum(Patterns(:,:,:,THETA_COMP,X_POL),3);
% Eph_sum_X = sum(Patterns(:,:,:,PHI_COMP,X_POL),3);
% 

%% plot 

% patterns
if polarization == 'L3'
    [EH,EV] = ThPh2L3(THETA*180/pi,PHI*180/pi,ETHETA_TOTAL,EPHI_TOTAL);
    norm=max(max([EH;EV]))
    plot3d(rad2deg(AZ),rad2deg(EL),db(EH/norm),'AZ','EL',[-40 0]);title('E_H total - Ludwig 3')
    plot3d(rad2deg(AZ),rad2deg(EL),db(EV/norm),'AZ','EL',[-40 0]);title('E_V total - Ludwig 3')
    errorEH = target_EL3co*N - abs(EH(TH_index,PH_index))
    errorEV = target_EL3x*N - abs(EV(TH_index,PH_index))
elseif polarization == 'L2'
    [EAz,EEl] = ThPh2L2I(THETA*180/pi,PHI*180/pi,ETHETA_TOTAL,EPHI_TOTAL);
    norm=max(max([EAz;EEl]))
    plot3d(rad2deg(AZ),rad2deg(EL),db(EAz/norm),'AZ','EL',[-40 0]);title('E_{Az} total - Ludwig 2I')
    plot3d(rad2deg(AZ),rad2deg(EL),db(EEl/norm),'AZ','EL',[-40 0]);title('E_{El} total - Ludwig 2I')
    errorEAz = target_EL2Az*N - abs(EAz(TH_index,PH_index))
    errorEEl = target_EL2El*N - abs(EEl(TH_index,PH_index))
end

%excitations
%pol 
normExc = max(max(ExcitationsCmplx));
figure
imagesc(reshape(db(ExcitationsCmplx(:,1)/normExc),NBR_ROWS,NBR_COLS));
colormap jet
colorbar
caxis([-30 0]);
ylabel('Rows')
xlabel('Cols')
title('X pol')

figure
imagesc(reshape(db(ExcitationsCmplx(:,2)/normExc),NBR_ROWS,NBR_COLS));
colormap jet
colorbar
caxis([-30 0]);
ylabel('Rows')
xlabel('Cols')
title('Y pol')

figure
bar(db(ExcitationsCmplx(:,:)/normExc),'stacked')
%plot(db(ExcitationsCmplx(:,:)/normExc))

% plot3d(rad2deg(THETA),rad2deg(PHI),db(ETHETA_TOTAL/norm),'Theta','Phi',[-40 0]);title('E_\theta total')
% plot3d(rad2deg(THETA),rad2deg(PHI),db(EPHI_TOTAL/norm),'Theta','Phi',[-40 0]);title('E_\phi total')
% plot3d(rad2deg(THETA),rad2deg(PHI),angle(ETHETA_TOTAL/norm)*180/pi,'Theta','Phi',[-180 180]);title('E_\theta total - phase')
% plot3d(rad2deg(THETA),rad2deg(PHI),angle(EPHI_TOTAL/norm)*180/pi,'Theta','Phi',[-180 180]);title('E_\phi total - phase')
% abs(ETHETA_TOTAL(TH_index,PH_index))
% abs(EPHI_TOTAL(TH_index,PH_index))
% eTH2plot = Patterns(:,:,4,THETA_COMP,X_POL)+Patterns(:,:,5,THETA_COMP,X_POL)+Patterns(:,:,6,THETA_COMP,X_POL);
% ePH2plot = Patterns(:,:,4,PHI_COMP,X_POL)+Patterns(:,:,5,PHI_COMP,X_POL)+Patterns(:,:,6,PHI_COMP,X_POL);
% norm=max(max(eTH2plot))
% plot3d(rad2deg(THETA),rad2deg(PHI),db(eTH2plot/norm),'Theta','Phi',[-40 0]);title('ETHETA')
% plot3d(rad2deg(THETA),rad2deg(PHI),db(ePH2plot/norm),'Theta','Phi',[-40 0]);title('EPHI')
% 
% plot3d(rad2deg(AZ),rad2deg(EL),10.*log10(abs(ETHETA_X).^2),'AZ','EL',[-40 40])
% plot3d(rad2deg(AZ),rad2deg(EL),10.*log10(abs(ETHETA_X_rot).^2),'AZ','EL',[-40 40])
% plot3d(rad2deg(AZ),rad2deg(EL),10.*log10(abs(ETHETA_Y).^2),'AZ','EL',[-40 40])
% plot3d(rad2deg(THETA),rad2deg(PHI),10.*log10(abs(boresight).^2),'AZ','EL',[-80 80])
