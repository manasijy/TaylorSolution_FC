function [gmatrix] = updatemat(spin,DCMtrx)
w12 = spin(1,2);
w31 = spin(3,1);
w23 = spin(2,3);


Net_rotaion = sqrt((w12*w12) + (w31*w31) + (w23*w23));
r1 = (-w23/Net_rotaion); r2 = (-w31/Net_rotaion); r3 = (-w12/Net_rotaion);
axis_angle_pair = [r1 r2 r3 Net_rotaion];
Rotation_Matrix = vrrotvec2mat(axis_angle_pair);%axang2rotm(axis_angle_pair);
New_rotation_matrix = Rotation_Matrix*DCMtrx';
% gmatrix = rotation('matrix',New_rotation_matrix); 

new_phi = real((radtodeg((acos(New_rotation_matrix(3,3))))));
new_phi1 = real((radtodeg((asin(New_rotation_matrix(3,1)/sin(degtorad(new_phi)))))));
new_phi2 = real((radtodeg((acos(New_rotation_matrix(2,3)/(sin(degtorad(new_phi))))))));
new_euler_angle = ([new_phi1 new_phi new_phi2]);
gmatrix = rotation('Euler',new_euler_angle*degree);


