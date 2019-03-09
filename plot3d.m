function [Plot1] = plot3d(X,Y,Data,xlab,ylab,colbar,visible)
if nargin < 7 
    visible=1;
end
%PLOT3D Produces a surface plot of "Data".
%   The function receives in the form of matrices X and Y variables of the 
%   mapped "Data" and produces a surface plot.
%
%   X           (m by n)        -   Y-axis coordinates
%   Y           (m by n)        -   Y-axis coordinates
%   Data        (m by n)        -   data to be plotted (z-axis)
%   xlab        (string vector) -   label for the x-axis
%   ylab        (string vector) -   label for the y-axis
%   colbar      (1 by 2)        -   data range for plotting
%
%   Example:
%
%   theta = linspace(-pi/2,pi/2,181);
%   phi = linspace(0,2*pi,181);
%   [THETA,PHI] = meshgrid(theta,phi);
%   u = sin(THETA).*cos(PHI);
%   v = sin(THETA).*sin(PHI);
%   w = cos(THETA);
%   
%   [Plot1] = plot3d(u,v,w,'u (rad)','v (rad)',[0 1]);
if visible
    Plot1 = figure('visible','on');
else
    Plot1 = figure('visible','off');
end
surf(X,Y,Data)
view([0 90])
shading interp
colormap jet
ylabel(ylab)
xlabel(xlab)
grid off
xlim([min(min(X)) max(max(X))])
ylim([min(min(Y)) max(max(Y))])
caxis(colbar)
set(gca,'fontsize',18)
colorbar;
zlim(colbar)

% saveas(Plot1,'Plot.fig')

end

