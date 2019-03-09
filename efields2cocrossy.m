function [u_co_2,u_cross_2,u_co_3,u_cross_3,u_co_2_dB,u_cross_2_dB,u_co_3_dB,u_cross_3_dB] = efields2cocrossy(THETA,PHI,Etheta,Ephi)

%efields2cocrossy Calculates Ludwig's Definition of Co- and 
%Cross-Polarization for an antenna that is polarized along the Y-axis.
%
%   THETA       (m by n)        -   Theta (radians) coordinates
%   PHI         (m by n)        -   Phi (radians) coordinates
%   Etheta      (m by n)        -   Electric Fields in Theta
%   EPhi        (m by n)        -   Electric Fields in Phi 
%
%   Example:
%
%   [u_co_2,u_cross_2,u_co_3,u_cross_3] = efields2cocrossx(THETA,PHI,Etheta,Ephi);

%%% Ludwig 2nd Definition of Radiation patterns 
u_co_2 = (cos(THETA).*sin(PHI).*Etheta+cos(PHI).*Ephi)./(sqrt(1-(sin(THETA).*sin(PHI)).^2));
u_cross_2 = (cos(PHI).*Etheta-cos(THETA).*sin(PHI).*Ephi)./(sqrt(1-(sin(THETA).*sin(PHI)).^2));

u_co_2_db = 10.*log10(abs(u_co_2).^2);
u_cross_2_db = 10.*log10(abs(u_cross_2).^2);

u_co_2_db(u_co_2_db==Inf) = 0;

Normalization_u_co_2 = max(max(u_co_2_db));

u_co_2_dB = u_co_2_db-Normalization_u_co_2;
u_cross_2_dB = u_cross_2_db-Normalization_u_co_2;

%%% Ludwig 3rd Definition of Radiation Patterns
u_co_3 = sin(PHI).*Etheta+cos(PHI).*Ephi;
u_cross_3 = cos(PHI).*Etheta-sin(PHI).*Ephi;

u_co_3_db = 10.*log10(abs(u_co_3).^2);
u_cross_3_db = 10.*log10(abs(u_cross_3).^2);

u_co_3_db(u_co_3_db==Inf) = 0;

Normalization_u_co_3 = max(max(u_co_3_db));

u_co_3_dB = u_co_3_db-Normalization_u_co_3;
u_cross_3_dB = u_cross_3_db-Normalization_u_co_3;

end

