function [U,V,W] = thph2uv(THETA,PHI)
%THPH2UV Calculates the Sperical Components UVW
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
%   [U,V,W] = thph2uv(THETA,PHI)
U = sin(THETA).*cos(PHI);
V = sin(THETA).*sin(PHI);
W = cos(THETA);

end

