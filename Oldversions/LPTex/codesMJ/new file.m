clear all;
clc;
desired_row = 0;
dummy1 = 0;
common_slip =[];
independent_slip_system = [];


%% Asking the euler angle input from User
prompt1 = 'Enter the Euler angles (in degrees) [phi1, phi, phi2] for example [35,20,20]:- ';
euler_angle = input(prompt1);
% prompt2 = 'Enter the strain increment (in decimal) for example [0.1]:- ';
% strain_increment = input(prompt2);
strain_increment =0.1;
% prompt2 = 'please specify the strain condition [exx exy exz;eyx eyy eyz; ezx ezy ezz] for example [1 0 0;0 0 0;0 0 -1];-';
% strain_external= input(prompt2);
strain_external = [1 0 0;0 0 0;0 0 -1];
%euler_angle

%% Calculating g_ij and transpose of it

%converting degree to radian and storing as phi1,phi,phi2
phi1 = degtorad(euler_angle(1));
phi = degtorad(euler_angle(2));
phi2 = degtorad(euler_angle(3));

g11 = ((cos(phi1)*cos(phi2))-sin(phi1)*sin(phi2)*cos(phi));
g12 = (sin(phi1)*cos(phi2))+(cos(phi1)*sin(phi2)*cos(phi));
g13 = sin(phi2)*sin(phi);

g21 = (-cos(phi1)*sin(phi2) - sin(phi1)*cos(phi2)*cos(phi));
g22 = (-sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(phi));
g23 = cos(phi2)*sin(phi);

g31 = sin(phi1)*sin(phi);
g32 = (-cos(phi1)*sin(phi));
g33 = cos(phi);

g = [g11 g12 g13;g21 g22 g23;g31 g32 g33];
z = [round(g13,1) round(g23,1) round(g33,1)];
x = [round(g11,1) round(g21,1) round(g31,1)];
y = [round(g12,1) round(g22,1) round(g32,1)];
g_inv=inv(g);
g_trans = g';
g_trans_inv=inv(g_trans);

%% Determinig the External Strain Matrix
strain_external;
%% calcualting the strains along crystallographic axis
s_in= g*strain_external*g_trans;

%% Table bishop-hill stress state

bh =    [1.00	-1.00	0.00	0.00	0.00	0.00
0.00	1.00	-1.00	0.00	0.00	0.00
-1.00	0.00	1.00	0.00	0.00	0.00
0.00	0.00	0.00	1.00	0.00	0.00
0.00	0.00	0.00	0.00	1.00	0.00
0.00	0.00	0.00	0.00	0.00	1.00
0.50	-1.00	0.50	0.00	0.50	0.00
0.50	-1.00	0.50	0.00	-0.50	0.00
-1.00	0.50	0.50	0.50	0.00	0.00
-1.00	0.50	0.50	-0.50	0.00	0.00
0.50	0.50	-1.00	0.00	0.00	0.50
0.50	0.50	-1.00	0.00	0.00	-0.50
0.50	0.00	-0.50	0.50	0.00	0.50
0.50	0.00	-0.50	-0.50	0.00	0.50
0.50	0.00	-0.50	0.50	0.00	-0.50
0.50	0.00	-0.50	-0.50	0.00	-0.50
0.00	-0.50	0.50	0.00	0.50	0.50
0.00	-0.50	0.50	0.00	-0.50	0.50
0.00	-0.50	0.50	0.00	0.50	-0.50
0.00	-0.50	0.50	0.00	-0.50	-0.50
-0.50	0.50	0.00	0.50	0.50	0.00
-0.50	0.50	0.00	-0.50	0.50	0.00
-0.50	0.50	0.00	0.50	-0.50	0.00
-0.50	0.50	0.00	-0.50	-0.50	0.00
0.00	0.00	0.00	0.50	0.50	-0.50
0.00	0.00	0.00	0.50	-0.50	0.50
0.00	0.00	0.00	-0.50	0.50	0.50
0.00	0.00	0.00	0.50	0.50	0.50
-1.00	1.00	0.00	0.00	0.00	0.00
0.00	-1.00	1.00	0.00	0.00	0.00
1.00	0.00	-1.00	0.00	0.00	0.00
0.00	0.00	0.00	-1.00	0.00	0.00
0.00	0.00	0.00	0.00	-1.00	0.00
0.00	0.00	0.00	0.00	0.00	-1.00
-0.50	1.00	-0.50	0.00	-0.50	0.00
-0.50	1.00	-0.50	0.00	0.50	0.00
1.00	-0.50	-0.50	-0.50	0.00	0.00
1.00	-0.50	-0.50	0.50	0.00	0.00
-0.50	-0.50	1.00	0.00	0.00	-0.50
-0.50	-0.50	1.00	0.00	0.00	0.50
-0.50	0.00	0.50	-0.50	0.00	-0.50
-0.50	0.00	0.50	0.50	0.00	-0.50
-0.50	0.00	0.50	-0.50	0.00	0.50
-0.50	0.00	0.50	0.50	0.00	0.50
0.00	0.50	-0.50	0.00	-0.50	-0.50
0.00	0.50	-0.50	0.00	0.50	-0.50
0.00	0.50	-0.50	0.00	-0.50	0.50
0.00	0.50	-0.50	0.00	0.50	0.50
0.50	-0.50	0.00	-0.50	-0.50	0.00
0.50	-0.50	0.00	0.50	-0.50	0.00
0.50	-0.50	0.00	-0.50	0.50	0.00
0.50	-0.50	0.00	0.50	0.50	0.00
0.00	0.00	0.00	-0.50	-0.50	0.50
0.00	0.00	0.00	-0.50	0.50	-0.50
0.00	0.00	0.00	0.50	-0.50	-0.50
0.00	0.00	0.00	-0.50	-0.50	-0.50];
%% Calculating taylor factor for single crystal

for i=1:1:56
    m(i) = round(sqrt(6)*( (bh(i,1)*s_in(2,2)) - ((bh(i,2)*s_in(1,1))) + ((bh(i,4)*2*s_in(2,3))) + ((bh(i,5)*2*s_in(3,1))) +((bh(i,6)*2*s_in(1,2)))),3);
    % 
end   
Taylor_factor = max(m);
for j=1:1:56
    if m(j)==max(m)
       desired_row = j;
       dummy1 = dummy1+1;

%% Defining Slip System

a = [1 1 1 0 1 -1; 1 1 1 -1 0 1; 1 1 1 1 -1 0; 1 1 1 0 -1 1; 1 1 1 1 0 -1; 1 1 1 -1 1 0];
b = [-1 -1 1 0 -1 -1; -1 -1 1 1 0 1; -1 -1 1 -1 1 0; -1 -1 1 0 1 1; -1 -1 1 -1 0 -1; -1 -1 1 1 -1 0];
c = [-1 1 1 0 1 -1; -1 1 1 1 0 1; -1 1 1 -1 -1 0; -1 1 1 0 -1 1; -1 1 1 -1 0 -1;-1 1 1 1 1 0];
d = [1 -1 1 0 -1 -1; 1 -1 1 -1 0 1; 1 -1 1 1 1 0; 1 -1 1 0 1 1; 1 -1 1 1 0 -1; 1 -1 1 -1 -1 0];

%% Defining operating slip system for all stress states

stress_state(:,:,1) = [a(1,:);a(5,:);b(1,:);b(5,:);c(1,:);c(5,:);d(1,:);d(5,:)];
stress_state(:,:,2) = [a(2,:);a(6,:);b(2,:);b(6,:);c(2,:);c(6,:);d(2,:);d(6,:)];
stress_state(:,:,3) = [a(4,:);a(3,:);b(4,:);b(3,:);c(4,:);c(3,:);d(4,:);d(3,:)];
stress_state(:,:,4) = [a(2,:);a(6,:);b(5,:);b(3,:);c(2,:);c(6,:);d(5,:);d(3,:)];
stress_state(:,:,5) = [a(4,:);a(3,:);b(1,:);b(6,:);c(1,:);c(6,:);d(4,:);d(3,:)];
stress_state(:,:,6) = [a(1,:);a(5,:);b(1,:);b(5,:);c(4,:);c(2,:);d(4,:);d(2,:)];
stress_state(:,:,7) = [a(5,:);a(3,:);b(1,:);b(5,:);c(1,:);c(5,:);d(5,:);d(3,:)];
stress_state(:,:,8) = [a(1,:);a(5,:);b(5,:);b(3,:);c(5,:);c(3,:);d(1,:);d(5,:)];
stress_state(:,:,9) = [a(4,:);a(2,:);b(4,:);b(3,:);c(4,:);c(2,:);d(4,:);d(3,:)];
stress_state(:,:,10) = [a(4,:);a(3,:);b(4,:);b(2,:);c(4,:);c(3,:);d(4,:);d(2,:)];
stress_state(:,:,11) = [a(1,:);a(6,:);b(1,:);b(6,:);c(2,:);c(6,:);d(2,:);d(6,:)];
stress_state(:,:,12) = [a(2,:);a(6,:);b(2,:);b(6,:);c(1,:);c(6,:);d(1,:);d(6,:)];
stress_state(:,:,13) = [a(1,:);a(6,:);b(1,:);b(5,:);c(2,:);c(6,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,14) = [a(1,:);a(5,:);b(1,:);b(6,:);d(2,:);d(6,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,15) = [a(2,:);a(6,:);c(1,:);c(6,:);d(1,:);d(5,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,16) = [b(2,:);b(6,:);c(1,:);c(5,:);d(1,:);d(5,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,17) = [a(5,:);a(3,:);b(1,:);b(5,:);d(4,:);d(3,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,18) = [a(1,:);a(5,:);b(5,:);b(3,:);c(4,:);c(3,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,19) = [a(4,:);a(3,:);c(1,:);c(5,:);d(5,:);d(3,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,20) = [b(4,:);b(3,:);c(5,:);c(3,:);d(1,:);d(5,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,21) = [a(4,:);a(2,:);c(2,:);c(6,:);d(4,:);d(3,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,22) = [a(4,:);a(3,:);b(2,:);b(6,:);d(4,:);d(2,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,23) = [a(2,:);a(6,:);b(4,:);b(3,:);c(4,:);c(2,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,24) = [b(4,:);b(2,:);c(4,:);c(3,:);d(2,:);d(6,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,25) = [a(4,:);a(2,:);c(1,:);c(6,:);d(5,:);d(3,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,26) = [a(1,:);a(6,:);b(5,:);b(3,:);c(4,:);c(2,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,27) = [a(5,:);a(3,:);b(1,:);b(6,:);d(4,:);d(2,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,28) = [b(1,:);b(5,:);c(2,:);c(6,:);d(4,:);d(3,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,29) = [a(4,:);a(2,:);b(4,:);b(2,:);c(4,:);c(2,:);d(4,:);d(2,:)];
stress_state(:,:,30) = [a(5,:);a(3,:);b(5,:);b(3,:);c(5,:);c(3,:);d(5,:);d(3,:)];
stress_state(:,:,31) = [a(1,:);a(6,:);b(1,:);b(6,:);c(1,:);c(6,:);d(1,:);d(6,:)];
stress_state(:,:,32) = [a(5,:);a(3,:);b(2,:);b(6,:);c(5,:);c(3,:);d(2,:);d(6,:)];
stress_state(:,:,33) = [a(1,:);a(6,:);b(4,:);b(3,:);c(4,:);c(3,:);d(1,:);d(6,:)];
stress_state(:,:,34) = [a(4,:);a(2,:);b(4,:);b(2,:);c(1,:);c(5,:);d(1,:);d(5,:)];
stress_state(:,:,35) = [a(2,:);a(6,:);b(4,:);b(2,:);c(4,:);c(2,:);d(2,:);d(6,:)];
stress_state(:,:,36) = [a(4,:);a(2,:);b(2,:);b(6,:);c(2,:);c(6,:);d(4,:);d(2,:)];
stress_state(:,:,37) = [a(1,:);a(5,:);b(1,:);b(6,:);c(1,:);c(5,:);d(1,:);d(6,:)];
stress_state(:,:,38) = [a(1,:);a(6,:);b(1,:);b(5,:);c(1,:);c(6,:);d(1,:);d(5,:)];
stress_state(:,:,39) = [a(4,:);a(3,:);b(4,:);b(3,:);c(5,:);c(3,:);d(5,:);d(3,:)];
stress_state(:,:,40) = [a(5,:);a(3,:);b(5,:);b(3,:);c(4,:);c(3,:);d(4,:);d(3,:)];
stress_state(:,:,41) = [a(4,:);a(3,:);b(4,:);b(2,:);c(5,:);c(3,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,42) = [a(4,:);a(2,:);b(4,:);b(3,:);d(5,:);d(3,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,43) = [a(5,:);a(3,:);c(4,:);c(3,:);d(4,:);d(2,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,44) = [b(5,:);b(3,:);c(4,:);c(2,:);d(4,:);d(2,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,45) = [a(2,:);a(6,:);b(4,:);b(2,:);d(1,:);d(6,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,46) = [a(4,:);a(2,:);b(2,:);b(6,:);c(1,:);c(6,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,47) = [a(1,:);a(6,:);c(4,:);c(2,:);d(2,:);d(6,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,48) = [b(1,:);b(6,:);c(2,:);c(6,:);d(4,:);d(2,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,49) = [a(1,:);a(5,:);c(5,:);c(3,:);d(1,:);d(6,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,50) = [a(1,:);a(6,:);b(5,:);b(3,:);d(1,:);d(5,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,51) = [a(5,:);a(3,:);b(1,:);b(6,:);c(1,:);c(5,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,52) = [b(1,:);b(5,:);c(1,:);c(6,:);d(5,:);d(3,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,53) = [a(1,:);a(5,:);c(4,:);c(3,:);d(2,:);d(6,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,54) = [a(4,:);a(3,:);b(2,:);b(6,:);c(1,:);c(5,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,55) = [a(2,:);a(6,:);b(4,:);b(3,:);d(1,:);d(5,:);0 0 0 0 0 0; 0 0 0 0 0 0];
stress_state(:,:,56) = [b(4,:);b(2,:);c(5,:);c(3,:);d(1,:);d(6,:);0 0 0 0 0 0; 0 0 0 0 0 0];

%% Coding to Output the Slip System 

possible_slip_system(:,:,dummy1) =stress_state(:,:,desired_row);
    end
end
slip_system_for_cal = possible_slip_system(:,:,1);
%% Finding the common solution for slip system 
% 
% if size(possible_slip_system,3)==2
% for l=1:1:8
%     for dummy2 = 1:1:8
%     if possible_slip_system(l,:,1)== possible_slip_system(dummy2,:,2)
%         common_slip = [common_slip;possible_slip_system(l,:,1)];
%     end
%     end
% end
% elseif size(possible_slip_system,3)==3
% for l=1:1:8
%     for dummy2 = 1:1:8
%         for dummy3 = 1:1:8
%     if possible_slip_system(l,:,1)== possible_slip_system(dummy2,:,2) & possible_slip_system(l,:,1)== possible_slip_system(dummy3,:,3)
%         common_slip = [common_slip;possible_slip_system(l,:,1)];
%     end
%         end
%     end
% end
% elseif size(possible_slip_system,3)==4
% for l=1:1:8
%     for dummy2 = 1:1:8
%         for dummy3 = 1:1:8
%             for dummy4 = 1:1:8
%     if possible_slip_system(l,:,1)== possible_slip_system(dummy2,:,2) & possible_slip_system(l,:,1)== possible_slip_system(dummy3,:,3) & possible_slip_system(l,:,1)== possible_slip_system(dummy4,:,4)
%         common_slip = [common_slip;possible_slip_system(l,:,1)];
%     end
%         end
%         end
%     end
% end
% prompt1 = 'Common Operative Slip System are:- ';
% common_slip;
% 
% 
% else
%     prompt1 = 'Possible Operative Slip System are:- ';
%     possible_slip_system;
%         
% end
% independent_slip_system = common_slip;
%% Calculating the shear for the case having unique solution



%% Generating all possible matrix
dummy_m = nchoosek(1:1:size(slip_system_for_cal,1),5);


for var1=1:1:size(dummy_m,1)
    generate_mat = [];
    for var2 = 1:1:size(dummy_m,2)
        generate_mat = [generate_mat;slip_system_for_cal(dummy_m(var1,var2),:)];
     
    end
    Three_DM(:,:,var1) = generate_mat;
    
end


%% calculating 5x5 matrix

for j=1:1:size(Three_DM,3)
  collect_column=[];  
for i=1:1:size(Three_DM,1)
column1 = (1/(sqrt(6)))*[Three_DM(i,5,j)*Three_DM(i,2,j),Three_DM(i,6,j)*Three_DM(i,3,j), Three_DM(i,5,j)*Three_DM(i,3,j)+Three_DM(i,6,j)*Three_DM(i,2,j), Three_DM(i,4,j)*Three_DM(i,3,j)+Three_DM(i,1,j)*Three_DM(i,6,j), Three_DM(i,4,j)*Three_DM(i,2,j)+Three_DM(i,1,j)*Three_DM(i,5,j)]';
collect_column = [collect_column column1];
end
fivebyfive_mat(:,:,j)=collect_column;
end

%% Calculating independent slip systems by eliminating depedent slip system
% D=0
% Shear = -ve
store_row_position = [];
collect_shear = [];
sd_shear = [];

for j=1:1:size(fivebyfive_mat,3)
    if det(fivebyfive_mat(:,:,j))==0
        continue
    end
    shear = inv(fivebyfive_mat(:,:,j))*[s_in(2,2);s_in(3,3);2*s_in(2,3);2*s_in(1,3);2*s_in(1,2)];
    if all(shear(:)>=0)
        
		if ismember(round(std(shear),4),round(sd_shear,4))
        continue
        end
	
		sd_shear = [sd_shear std(shear)];
        collect_shear=[collect_shear shear];
        store_row_position = [store_row_position j];
        
    else
        continue
    end
end

%% Finding the Minimum Standard deviation for Shear

[min_sd index_min_sd] = min(sd_shear);
   
location_of_ss = store_row_position(index_min_sd);

independent_slip_system = Three_DM(:,:,location_of_ss);


%% Forming 3x3 matrix by calculating diadic prodcut

for i=1:1:size(independent_slip_system,1)
three_by_three(:,:,i) = [independent_slip_system(i,4);independent_slip_system(i,5); independent_slip_system(i,6)]*[independent_slip_system(i,1) independent_slip_system(i,2) independent_slip_system(i,3)];
end

%% finding e12, e21, e13, e31, e23 and e32
e12=0; e21=0; e13=0; e31=0; e23=0; e32 = 0;

% finding e12
for i = 1:1:size(three_by_three,3)
dummy2 = three_by_three(1,2,i)*collect_shear(i)*strain_increment;
e12 = dummy2+e12;
end


% finding e21
for i = 1:1:size(three_by_three,3)
dummy2 = three_by_three(2,1,i)*collect_shear(i)*strain_increment;
e21 = dummy2+e21;
end


% finding e13
for i = 1:1:size(three_by_three,3)
dummy2 = three_by_three(1,3,i)*collect_shear(i)*strain_increment;
e13 = dummy2+e13;
end


% finding e31
for i = 1:1:size(three_by_three,3)
dummy2 = three_by_three(3,1,i)*collect_shear(i)*strain_increment;
e31 = dummy2+e31;
end


% finding e23
for i = 1:1:size(three_by_three,3)
dummy2 = three_by_three(2,3,i)*collect_shear(i)*strain_increment;
e23 = dummy2+e23;
end


% finding e32
for i = 1:1:size(three_by_three,3)
dummy2 = three_by_three(3,2,i)*collect_shear(i)*strain_increment;
e32 = dummy2+e32;
end


%% calculating W12, W31 and W23

w12 = (e12-e21)/4.88;
w13 = (e13-e31)/4.88;
w23 = (e23-e32)/4.88;
w21 = (e21-e12)/4.88;
w31 = (e31-e13)/4.88;
w32 = (e32-e23)/4.88;
spin= [0 w12 w13;w21 0 w23;w31 w32 0];
quaternion= matrix2quat(spin);
rotation_matrix=quat2rotm(quaternion);
New_rotation_matrix= rotation_matrix*g;

% Net_rotaion = sqrt((w12*w12) + (w31*w31) + (w23*w23));
% r1 = (w23/Net_rotaion); r2 = (w31/Net_rotaion); r3 = (w12/Net_rotaion);
% axis_angle_pair = [r1 r2 r3 Net_rotaion];
% Rotation_Matrix = axang2rotm(axis_angle_pair);
% New_rotation_matrix = Rotation_Matrix*g;

% new_phi = radtodeg(acos(New_rotation_matrix(3,3)));
% new_phi1 = radtodeg(asin(New_rotation_matrix(3,1)/sin(degtorad(new_phi))));
% new_phi2 = radtodeg(acos(New_rotation_matrix(2,3)/sin(degtorad(new_phi))));
% new_euler_angle = [new_phi1 new_phi new_phi2];





