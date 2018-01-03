function [g] = euler2gmatrix(EU)

% EU is [phi1,phi,phi2]. It is the set of bunge euler angles and it is
% input the function. Following code generates the orientation matrix g for
% this orientation
% psi theta and phi are the kocks sphereical angles and these are related
% to bunge euler angles by the following equations.


psi = EU{1,1}-pi/2;
theta = EU{1,2};
phi = pi/2 - EU{1,3};

g11 = -sind(phi).*sind(psi)-cosd(phi).*cosd(psi).*cosd(theta);
g12 = sind(phi).*cosd(psi)-cosd(phi).*sind(psi).*cosd(theta);
g13 = cosd(phi).*sind(theta);
g21 = cosd(phi).*sind(psi)-sind(phi).*cosd(psi).*cosd(theta);
g22 = -cosd(phi).*cosd(psi)-sind(phi).*sind(psi).*cosd(theta);
g23 = sind(phi).*sind(theta);
g31 = cosd(psi).*sind(theta);
g32 = sind(psi).*sind(theta);
g33 = cosd(theta);
g = [g11,g12,g13;g21,g22,g23;g31,g32,g33];