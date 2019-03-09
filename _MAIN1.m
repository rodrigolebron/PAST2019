close all
clear all
clc

theta = -45*pi/180;
phi = 45*pi/180;

[AZ,EL] = thph2azel(theta,phi);

AZ*180/pi
EL*180/pi

[THETA,PHI] = azel2thph(AZ,EL);

THETA*180/pi
PHI*180/pi