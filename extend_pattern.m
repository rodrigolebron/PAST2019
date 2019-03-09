function [THETA_NEW,PHI_NEW,PATTERN_NEW] = extend_pattern(THETA,PHI,PATTERN,k,Figures)
%extend_pattern Summary of this function goes here
%   Detailed explanation goes here

if THETA(1,1) ~= THETA(2,1)
    
    THETA = THETA';
    PHI = PHI';
    PATTERN = PATTERN';
    disp('Data was re-arranged.');
    
else
    
    disp('Data was NOT re-arranged.');
    
end
    
%%% Interpolating Pattern

THETA_NEW = interp2(THETA,k);
PHI_NEW = interp2(PHI,k);

PATTERN_NEW = interp2(THETA,PHI,PATTERN,THETA_NEW,PHI_NEW);

%%% Plots

if Figures == 1

    figure
    imagesc(abs(PATTERN))
    colormap jet
    
    figure
    imagesc(abs(PATTERN_NEW))
    colormap jet

else

    disp('Plots Disabled')
    
end

end