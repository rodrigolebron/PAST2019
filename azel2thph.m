function [THETA,PHI] = azel2thph(AZ,EL)

%AZEL2THPH Calculates Coordinate System THETA/PHI
%
%   AZ         (m by n)
%   EL           (m by n)
%
%   Example:
%
%   az = linspace(-pi/2,pi/2,181);
%   el = linspace(-pi/2,pi/2,181);
%   [AZ,EL] = meshgrid(az,el);
% 
%   [THETA,PHI] = azel2thph(AZ,EL)

THETA = acos(cos(EL).*cos(AZ));
PHI = atan2(tan(EL),sin(AZ));

end

