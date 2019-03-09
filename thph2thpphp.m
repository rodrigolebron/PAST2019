function [THETA_P,PHI_P,Htheta_P,Hphi_P] = thph2thpphp(alpha,beta,gamma,THETA,PHI,Ephi,Etheta,Figures)
%THPH2THPPHP Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preamble
% Etheta_Ephi = csvread('Etheta_Ephi_Xpol.csv',1,0);

% THETA = deg2rad(reshape(Etheta_Ephi(:,1),181,[]));
% PHI = deg2rad(reshape(Etheta_Ephi(:,2),181,[]));

M = length(THETA);
N = length(PHI);

% Ephi = (reshape(Etheta_Ephi(:,3),181,[])+1i.*reshape(Etheta_Ephi(:,4),181,[]))./1000;
% Etheta = (reshape(Etheta_Ephi(:,5),181,[])+1i.*reshape(Etheta_Ephi(:,6),181,[]))./1000;

H = {zeros(M,N); Etheta; Ephi};

f1 = 0;
f2 = 0;
f3 = 0;

r = 1;
%r_p = r-f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Building the "Rotation" Matrix

A11 = ones(M,N).*(cos(gamma).*cos(alpha)-sin(gamma).*cos(beta).*sin(alpha));
A12 = ones(M,N).*(cos(gamma).*sin(alpha)+sin(gamma).*cos(beta).*cos(alpha));
A13 = ones(M,N).*(sin(gamma).*sin(beta));

A21 = ones(M,N).*(-sin(gamma).*cos(alpha)-cos(gamma).*cos(beta).*sin(alpha));
A22 = ones(M,N).*(-sin(gamma).*sin(alpha)+cos(gamma).*cos(beta).*cos(alpha));
A23 = ones(M,N).*(cos(gamma).*sin(beta));

A31 = ones(M,N).*(sin(beta).*sin(alpha));
A32 = ones(M,N).*(-sin(beta).*cos(alpha));
A33 = ones(M,N).*(cos(beta));

A = {A11 A12 A13 ; A21 A22 A23 ; A31 A32 A33};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculating "Theta' and Phi' Coordinate System"
THP_PHP1 = zeros(M,N);
THP_PHP2 = zeros(M,N);
THP_PHP3 = zeros(M,N);

TH_PH = {r.*sin(THETA).*cos(PHI)-f1, zeros(M,N), zeros(M,N) ; r.*sin(THETA).*sin(PHI)-f2, zeros(M,N), zeros(M,N) ; r.*cos(THETA)-f3, zeros(M,N), zeros(M,N)};


for m = 1 : M
    for n = 1 : N
        
        A_Matrix = cellfun(@(c)c(m,n), A);
        TH_PH_Matrix = cellfun(@(c)c(m,n), TH_PH);
        
        Product = A_Matrix*TH_PH_Matrix;
        
        THP_PHP1(m,n) = Product(1);
        THP_PHP2(m,n) = Product(2);
        THP_PHP3(m,n) = Product(3);
        
    end
end

THETA_P = acos(THP_PHP3);
PHI_P = atan2(THP_PHP2,THP_PHP1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Building "XYZ to Spherical Transform" Matrix
B11 = sin(THETA).*cos(PHI);
B12 = sin(THETA).*sin(PHI);
B13 = cos(THETA);

B21 = cos(THETA).*cos(PHI);
B22 = cos(THETA).*sin(PHI);
B23 = -sin(THETA);

B31 = -sin(PHI);
B32 = cos(PHI);
B33 = zeros(M,N);

B = {B11 B12 B13 ; B21 B22 B23 ; B31 B32 B33};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Building "Spherical to XYZ Transform" Matrix

C11 = sin(THETA_P).*cos(PHI_P);
C12 = cos(THETA_P).*cos(PHI_P);
C13 = -sin(PHI_P);

C21 = sin(THETA_P).*sin(PHI_P);
C22 = cos(THETA_P).*sin(PHI_P);
C23 = cos(PHI_P);

C31 = cos(THETA_P);
C32 = -sin(THETA_P);
C33 = zeros(M,N);

C = {C11 C12 C13 ; C21 C22 C23 ; C31 C32 C33};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculating H-r', H-theta', and H-phi'

Hr_P = zeros(M,N);
Htheta_P = zeros(M,N);
Hphi_P = zeros(M,N);


for m = 1 : M
    for n = 1 : N
        
        B_Matrix = cellfun(@(c)c(m,n), B);
        A_Matrix = cellfun(@(c)c(m,n), A');
        C_Matrix = cellfun(@(c)c(m,n), C);
        H_Matrix = cellfun(@(c)c(m,n), H);
        
        Product = B_Matrix*A_Matrix*C_Matrix;
        
        Product_Inv = inv(Product);
        
        H_P  = Product_Inv\H_Matrix;
        
        Hr_P(m,n) = H_P(1);
        Htheta_P(m,n) = H_P(2);
        Hphi_P(m,n) = H_P(3);
        
    end
end

if Figures == 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting Un-Rotated Version

Etotal = sqrt(abs(Etheta).^2+abs(Ephi).^2);

u = Etotal.*sin(THETA).*cos(PHI);
v = Etotal.*sin(THETA).*sin(PHI);
w = Etotal.*cos(THETA);

figure
surf(u,v,w,'FaceColor','interp','FaceLighting','phong')
camlight right
xlabel('X')
ylabel('Y')
zlabel('Z')
shading interp
axis square

xlim([-15 15])
ylim([-15 15])
zlim([-15 15])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting Rotated Version

Etotal_P = sqrt(abs(Htheta_P).^2+abs(Hphi_P).^2);

u = Etotal_P.*sin(THETA_P).*cos(PHI_P);
v = Etotal_P.*sin(THETA_P).*sin(PHI_P);
w = Etotal_P.*cos(THETA_P);

figure
surf(u,v,w,'FaceColor','interp','FaceLighting','phong')
camlight right
xlabel('X')
ylabel('Y')
zlabel('Z')
shading interp
axis square

xlim([-15 15])
ylim([-15 15])
zlim([-15 15])

else
    
    disp('Plots Disabled')
    
end

end

