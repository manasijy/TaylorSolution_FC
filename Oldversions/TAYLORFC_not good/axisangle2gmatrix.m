function [g] = axisangle2gmatrix(rv)
%rv is rotation_vector
n1 = rv(1);
n2 = rv(2);
n3 = rv(3);
omega = rv(4); % in radians

g(1,1) = (1-n1^2)*cos(omega) +n1^2;
g(1,2) = n1*n2*(1-cos(omega))+n3*sin(omega);
g(1,3) = n1*n3*(1-cos(omega)) -n2*sin(omega);
g(2,1) = n2*n1*(1-cos(omega))-n3*sin(omega);
g(2,2) = (1-n2^2)*cos(omega) +n2^2;
g(2,3) = n2*n3*(1-cos(omega))-n1*sin(omega);
g(3,1) = n3*n1*(1-cos(omega))+n2*sin(omega);
g(3,2) = n3*n2*(1-cos(omega))-n1*sin(omega);
g(3,3) = (1-n3^2)*cos(omega) +n3^2;