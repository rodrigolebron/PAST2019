function [AZ,EL] = thph2azel(THETA,PHI)

%THPH2AZEL Calculates Coordinate System AZ/EL
%
%   THETA         (m by n)
%   PHI           (m by n)
%
%   Example:
%
%   theta = linspace(-pi/2,pi/2,181);
%   phi = linspace(0,2*pi,181);
%   [THETA,PHI] = meshgrid(theta,phi);
%
%   [AZ,EL] = thph2azel(THETA,PHI)

EL = asin(sin(PHI).*sin(THETA));
AZ = atan(cos(PHI).*tan(THETA));

end

