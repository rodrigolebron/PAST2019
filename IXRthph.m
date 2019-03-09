function [IXR_ThPh_dB] = IXRthph(THETA,PHI,ETHETA_X,EPHI_X,ETHETA_Y,EPHI_Y,Figures)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[AZ,EL] = thph2azel(THETA,PHI);

Condition_Number_Matrix_ThPh = zeros(181,181);

for m = 1 : 181
    for n = 1 : 181

        Jones_Matrix_ThPh = [EPHI_X(m,n),EPHI_Y(m,n);ETHETA_X(m,n),ETHETA_Y(m,n)];
        Condition_Number_Matrix_ThPh(m,n) = cond(Jones_Matrix_ThPh);         
        
    end
end

IXR_ThPh = ((Condition_Number_Matrix_ThPh+1)./(Condition_Number_Matrix_ThPh-1)).^2;

IXR_ThPh_dB = 10.*log10(IXR_ThPh);

if Figures == 1

    plot3d(rad2deg(AZ),rad2deg(EL),IXR_ThPh_dB,'Azimuth','Elevation',[0 60]);

else
    
    disp('Plots Disabled')
    
end

end

