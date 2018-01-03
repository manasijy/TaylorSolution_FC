function [g] = euler2gmatrix(EU)

% EU is [phi1,phi,phi2]. It is the set of bunge euler angles and it is
% input the function. Following code generates the orientation matrix g for
% this orientation
% psi theta and phi are the kocks sphereical angles and these are related
% to bunge euler angles by the following equations.

psi = EU(1)-pi/2;
theta = EU(2);
phi = pi/2 - EU(3);

g(1,1) = -sind(phi)*sind(psi)-cosd(phi)*cosd(psi)*cosd(theta);
g(1,2) = sind(phi)*cosd(psi)-cosd(phi)*sind(psi)*cosd(theta);
g(1,3) = cosd(phi)*sind(theta);
g(2,1) = cosd(phi)*sind(psi)-sind(phi)*cosd(psi)*cosd(theta);
g(2,2) = -cosd(phi)*cosd(psi)-sind(phi)*sind(psi)*cosd(theta);
g(2,3) = sind(phi)*sind(theta);
g(3,1) = cosd(psi)*sind(theta);
g(3,2) = sind(psi)*sind(theta);
g(3,3) = cosd(theta);