function [Array_Pattern,Array_Pattern_dB] = A2DS(Frequency,Nx,Ny,Dx,Dy,Theta_Scan,Phi_Scan,THETA,PHI)

%   Example: A2DS(3e9,8,8,50e-3,50e-3,30,30,THETA_P,PHI_P)

Resolution = size(THETA);   
 
c = 3e8;                        %Speed of Light [meters/seconds]
 
Theta_Scan = Theta_Scan*pi/180; %Scan Angle in Theta [radians]
Phi_Scan = Phi_Scan*pi/180;     %Scan Angle in Phi [radians]

THETA = THETA*pi/180;           %Theta Grid [degrees]
PHI = PHI*pi/180;               %Theta Grid [degrees]
     
lambda = c/Frequency;           %Wavelength [meters]
k = 2*pi/lambda;                %Wave number [radians/meters]
 
bx = -k.*Dx.*sin(Theta_Scan).*cos(Phi_Scan);     %Propagation Constant in the X axis
by = -k.*Dy.*sin(Theta_Scan).*sin(Phi_Scan);     %Propagation Constant in the Y axis
 
Array_Pattern = ones(Resolution);
 
Pattern_Add_X = 0;  %Pattern addition for X-Plane
Pattern_Add_Y = 0;  %Pattern addition for Y-Plane
 
for m = 1 : Nx
    Pattern_Add_X = Pattern_Add_X + Array_Pattern.*exp(1i.*(m-1).*(k.*Dx.*sin(THETA).*cos(PHI)+bx));
end
 
for n = 1 : Ny
    Pattern_Add_Y = Pattern_Add_Y + Array_Pattern.*exp(1i.*(n-1).*(k.*Dy.*sin(THETA).*sin(PHI)+by));
end
  
Array_Pattern = Pattern_Add_X.*Pattern_Add_Y;

Array_Pattern_dB = 10.*log10(abs(Array_Pattern).^2)-max(max(10.*log10(abs(Array_Pattern).^2)));
 
end