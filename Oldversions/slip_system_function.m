%% Slip system selection

% This program takes tensile direction in xtal reference frame as input and calculates schmid factor for all slip systems
% Then it selects slip systems with the highest schmidt factor

function [favourableSS, n_favourableSS, m_max] = slip_system_function(d)

% fprintf('Please enter the loading direction \n')
% O= input('');         % crystal_orientation array
% d = sscanf(O,'%d');
% d = d.' ;
m = zeros(1,24);
mod_n=1.7321;    
mod_b=1.4142;

% SS means slip system.

SlipSystem = struct('n',zeros(1,3),'b',zeros(1,3));

SlipSystem(1).n=[1,1,1]; SlipSystem(1).b =[0,-1,1];
SlipSystem(2).n=[1,1,1]; SlipSystem(2).b =[-1,0,1];
SlipSystem(3).n=[1,1,1]; SlipSystem(3).b =[-1,1,0];
SlipSystem(4).n=[-1,1,1]; SlipSystem(4).b =[1,1,0];
SlipSystem(5).n=[-1,1,1];SlipSystem(5).b =[1,0,1];
SlipSystem(6).n=[-1,1,1];SlipSystem(6).b =[0,-1,1];
SlipSystem(7).n=[1,-1,1]; SlipSystem(7).b =[1,1,0];
SlipSystem(8).n=[1,-1,1]; SlipSystem(8).b =[-1,0,1];
SlipSystem(9).n=[1,-1,1]; SlipSystem(9).b =[0,1,1];
SlipSystem(10).n=[-1,-1,1];SlipSystem(10).b =[-1,1,0];
SlipSystem(11).n=[-1,-1,1]; SlipSystem(11).b=[1,0,1];
SlipSystem(12).n=[-1,-1,1]; SlipSystem(12).b =[0,1,1];

SlipSystem(13).n=[1,1,1]; SlipSystem(13).b =[0,1,-1];
SlipSystem(14).n=[1,1,1]; SlipSystem(14).b =[1,0,-1];
SlipSystem(15).n=[1,1,1]; SlipSystem(15).b =[1,-1,0];
SlipSystem(16).n=[-1,1,1]; SlipSystem(16).b =[-1,-1,0];
SlipSystem(17).n=[-1,1,1];SlipSystem(17).b =[-1,0,-1];
SlipSystem(18).n=[-1,1,1];SlipSystem(18).b =[0,1,-1];
SlipSystem(19).n=[1,-1,1]; SlipSystem(19).b =[-1,-1,0];
SlipSystem(20).n=[1,-1,1]; SlipSystem(20).b =[1,0,-1];
SlipSystem(21).n=[1,-1,1]; SlipSystem(21).b =[0,-1,-1];
SlipSystem(22).n=[-1,-1,1];SlipSystem(22).b =[1,-1,0];
SlipSystem(23).n=[-1,-1,1]; SlipSystem(23).b=[-1,0,-1];
SlipSystem(24).n=[-1,-1,1]; SlipSystem(24).b =[0,-1,-1];

for i=1:1:24 
    
    n_vector = SlipSystem(i).n;
    b_vector = SlipSystem(i).b;
    
    mod_d=sqrt(d(1)^2+ d(2)^2+ d(3)^2);
    cos_lambeda= dot(d,b_vector)/(mod_b*mod_d);
    cos_phi= dot(d,n_vector)/(mod_n*mod_d);
    m(i) = cos_lambeda*cos_phi;
end

[m_max, ss_number] = max(m);

% m_max=0;
j = 1;
for i=1:1:24
    if m(i) == m_max
        max_position(j) = i;
        j = j+1;
%         m_max = m(i);
%         ss_number=i;
    end
end
n_favourableSS = length(max_position);
for k =1:1:n_favourableSS
    favourableSS(k) = SlipSystem(max_position(k));
%     favourableSS(k) = SlipSystem(k);
    fprintf('\n SS %8d B: %8d %d %d \t N: %d %d %d \n',max_position(k),SlipSystem(max_position(k)).b,SlipSystem(max_position(k)).n);
end

fprintf('Schmidt factor m = %f \n',m_max);


end
    
