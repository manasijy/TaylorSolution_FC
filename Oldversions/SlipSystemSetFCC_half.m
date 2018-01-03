%% This program creates slip system set including schimd factor for each system
% Both positive and negative SS are taken in to consideration

CS = crystalSymmetry('-43m');

SlipSystem = struct('n',zeros(1,3),'b',zeros(1,3));
SlipSystem(1).n=[1,1,1]; SlipSystem(1).b =[0,-1,1];%
SlipSystem(2).n=[1,1,1]; SlipSystem(2).b =[-1,0,1];%a2
SlipSystem(3).n=[1,1,1]; SlipSystem(3).b =[-1,1,0];
SlipSystem(4).n=[-1,1,1]; SlipSystem(4).b =[1,1,0];
SlipSystem(5).n=[-1,1,1];SlipSystem(5).b =[1,0,1];
SlipSystem(6).n=[-1,1,1];SlipSystem(6).b =[0,-1,1];
SlipSystem(7).n=[1,-1,1]; SlipSystem(7).b =[1,1,0];
SlipSystem(8).n=[1,-1,1]; SlipSystem(8).b =[-1,0,1];
SlipSystem(9).n=[1,-1,1]; SlipSystem(9).b =[0,1,1];
SlipSystem(10).n=[1,1,-1];SlipSystem(10).b =[-1,1,0];
SlipSystem(11).n=[1,1,-1]; SlipSystem(11).b=[1,0,1];
SlipSystem(12).n=[1,1,-1]; SlipSystem(12).b =[0,1,1];%b1


% SlipSystem(13).n=[1,1,1]; SlipSystem(13).b =[0,1,-1];%a1
% SlipSystem(14).n=[1,1,1]; SlipSystem(14).b =[1,0,-1];%
% SlipSystem(15).n=[1,1,1]; SlipSystem(15).b =[1,-1,0];%a3
% SlipSystem(16).n=[-1,1,1]; SlipSystem(16).b =[-1,-1,0];
% SlipSystem(17).n=[-1,1,1];SlipSystem(17).b =[-1,0,-1];
% SlipSystem(18).n=[-1,1,1];SlipSystem(18).b =[0,1,-1];
% SlipSystem(19).n=[1,-1,1]; SlipSystem(19).b =[-1,-1,0];
% SlipSystem(20).n=[1,-1,1]; SlipSystem(20).b =[1,0,-1];
% SlipSystem(21).n=[1,-1,1]; SlipSystem(21).b =[0,-1,-1];
% SlipSystem(22).n=[1,1,-1];SlipSystem(22).b =[1,-1,0];%b3
% SlipSystem(23).n=[1,1,-1]; SlipSystem(23).b=[-1,0,-1];%b2
% SlipSystem(24).n=[1,1,-1]; SlipSystem(24).b =[0,-1,-1];

%% NonOctahedral Slip Systems: (001)[110]

% % SlipSystem(13).n=[1,0,0]; SlipSystem(25).b =[0,1,1];
% % SlipSystem(14).n=[1,0,0]; SlipSystem(26).b =[0,-1,1];
% % SlipSystem(15).n=[0,1,0]; SlipSystem(27).b =[1,0,1];
% % SlipSystem(16).n=[0,1,0]; SlipSystem(28).b =[-1,0,1];
% % SlipSystem(17).n=[0,0,1];SlipSystem(29).b =[1,1,0];
% % SlipSystem(18).n=[0,0,1];SlipSystem(30).b =[-1,1,0];

% SlipSystem(25).n=[1,0,0]; SlipSystem(25).b =[0,1,1];
% SlipSystem(26).n=[1,0,0]; SlipSystem(26).b =[0,-1,1];
% SlipSystem(27).n=[0,1,0]; SlipSystem(27).b =[1,0,1];
% SlipSystem(28).n=[0,1,0]; SlipSystem(28).b =[-1,0,1];
% SlipSystem(29).n=[0,0,1];SlipSystem(29).b =[1,1,0];
% SlipSystem(30).n=[0,0,1];SlipSystem(30).b =[-1,1,0];

% SlipSystem(31).n=[1,0,0]; SlipSystem(31).b =[0,-1,-1];
% SlipSystem(32).n=[1,0,0]; SlipSystem(32).b =[0,1,-1];
% SlipSystem(33).n=[0,1,0]; SlipSystem(33).b =[-1,0,-1];
% SlipSystem(34).n=[0,1,0]; SlipSystem(34).b =[1,0,-1];
% SlipSystem(35).n=[0,0,1];SlipSystem(35).b =[-1,-1,0];
% SlipSystem(36).n=[0,0,1];SlipSystem(36).b =[1,-1,0];

for i = 1:1:12;%36

n = SlipSystem(i).n;
b = SlipSystem(i).b;
hkl = Miller(n(1),n(2),n(3),CS);
uvw = Miller(b(1),b(2),b(3),CS);

R = SchmidTensor(hkl,uvw);  %,'generalized');
% Rmatrix = matrix(R);
SlipSystem(i).SchmidT = R; %.M;

end
SlipSystemSet_FCC_half = SlipSystem;
save('SS_FCC_half.mat','SlipSystemSet_FCC_half');
% Rmatrix = matrix(R);
% m = vector3d.SchmidTensor(n,b);
% matrix(R)
% tau = double(EinsteinSum(R,[-1,-2],sigma,[-1,-2],'name','Schmid factor'));