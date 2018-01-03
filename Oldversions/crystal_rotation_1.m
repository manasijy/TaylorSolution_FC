%% Grain rotation 
% This program takes tensile direction in xtal reference frame and slip
% systems with corresponding shears(d_gamma) as input and then plots the  
% axis rotation wrt crystal   


B = zeros(1,3);
N = zeros(1,3);
fprintf('Please enter the loading direction e.g. ''1 1 1'' \n')
O= input('');         % crystal_orientation array
d0 = sscanf(O,'%d'); % d becomes column vector

fprintf('Please enter the slip plane indices \n')
On= input('');         % crystal_orientation array
SlipSystem(1).n = sscanf(On,'%d'); % d becomes column vector

fprintf('Please enter the slip direction indices \n')
Ob= input('');         % crystal_orientation array
SlipSystem(1).b = sscanf(Ob,'%d'); 

% fprintf('Please enter the d_gamma for this slip system \n')
% Og= input('');         % crystal_orientation array
% d_gamma = sscanf(Og,'%d'); 


CS = crystalSymmetry('-43m');
D = Miller(0,0,1,CS);
    

    N = SlipSystem(1).n;
    N = N/norm(N);
    B = SlipSystem(1).b;
    B = B/norm(B);

        d = d0;
        for j=1:1:21
            gamma=0.01*(j-1);
            for i=1:1:3
            d(i)= d(i)+ gamma*(dot(d,N))*B(i);
            end

         D(j) = Miller(d(1),d(2),d(3),CS);
        end      

h = plotstereogram;
plot_plasticSpin(h,D);
