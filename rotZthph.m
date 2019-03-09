function [ETHETA_X_rot,EPHI_X_rot,ETHETA_Y_rot,EPHI_Y_rot] = rotZthph(rot_angle,ETHETA_X,EPHI_X,ETHETA_Y,EPHI_Y)
%assumes PHI is changin in cols e.g. PHI=[0,1,2,3,...;0,1,2,3,...]
%and 1 deg spacing
%spherical chickes in vacuum rotation
%rotates around Z-axis an angle alpha

col_shift = floor(rot_angle);

ETHETA_X_rot = circshift(ETHETA_X,col_shift,2); 
EPHI_X_rot = circshift(EPHI_X,col_shift,2); 
ETHETA_Y_rot = circshift(ETHETA_Y,col_shift,2); 
EPHI_Y_rot = circshift(EPHI_Y,col_shift,2); 

if col_shift > 0
    
    ETHETA_X_rot(:,1:col_shift) = flipud(ETHETA_X_rot(:,1:col_shift).*exp(-1j*pi));
    EPHI_X_rot(:,1:col_shift) = flipud(EPHI_X_rot(:,1:col_shift).*exp(-1j*pi));
    ETHETA_Y_rot(:,1:col_shift) = flipud(ETHETA_Y_rot(:,1:col_shift).*exp(-1j*pi));
    EPHI_Y_rot(:,1:col_shift) = flipud(EPHI_Y_rot(:,1:col_shift).*exp(-1j*pi));
    
    else if col_shift < 0
 
        ETHETA_X_rot(:,length(ETHETA_X)-abs(col_shift+1):end) = flipud(ETHETA_X_rot(:,length(ETHETA_X)-abs(col_shift+1):end).*exp(-1j*pi));
        EPHI_X_rot(:,length(ETHETA_X)-abs(col_shift+1):end) = flipud(EPHI_X_rot(:,length(ETHETA_X)-abs(col_shift+1):end).*exp(-1j*pi));
        ETHETA_Y_rot(:,length(ETHETA_X)-abs(col_shift+1):end) = flipud(ETHETA_Y_rot(:,length(ETHETA_X)-abs(col_shift+1):end).*exp(-1j*pi));
        EPHI_Y_rot(:,length(ETHETA_X)-abs(col_shift+1):end) = flipud(EPHI_Y_rot(:,length(ETHETA_X)-abs(col_shift+1):end).*exp(-1j*pi));
    
        else
        
            disp('No rotation occurred.')
        
    end
end
    
    
end