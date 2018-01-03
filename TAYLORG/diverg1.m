
clear;
%% 
% Here we are creating a grid for a constant phi2.

v = 0:5:25;
[X, Y] = meshgrid(v);
%% 
% Now converting the grid to phi1,phi coordinate arrays.

phi1 = reshape(X,[],1);
phi  = reshape(Y,[],1);
%% 
% Consider phi2 = 45 section

phi2 = 45*ones(length(phi1),1);
p0 = [phi1, phi, phi2];
%% 
% Now calculate the orientation matrix

gmatrix0 = rotation('Euler',p0*degree);
% gmatrix0 = rotation('Euler',[phi1 , phi , phi2]*degree);
%% 
% f (full), p (only partial), fp (full+partial), ft(full+twin)
%%

slip_type = 'f';
%% 
%  0 (Full), 1 (Lath), 2 (PanCake)

constraints = 0;
max_strain = .02;
incr_e = 0.02;
%% 
% 0 (minVar), 1(maxVar), 2(minPlasticSpin), 3(wintenberger) 


criteria = 0;

tic
%%
%% Initializing the cost matrix according to type of slip mode
c = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
% full slip: n = 1-12,25-36 
% partial slip: n=13-24,37-48
% twin:n=13-24

switch slip_type
    case {'F', 'f'} % only perfect slip
        n = [13:24,37:48];
        c(n)=1e4; % All partial modes are made costlier        
    case {'P','p'} % only partial slip
        n = [1:12,25:36];
        c(n)=1e4; % All perfect modes are made costlier
    case {'FP','fp','Fp','fP'} %both perfect and partial slips
        
    case {'FT','ft','Ft','fT'}
        n = 37:48;
        c(n)=1e4; % All - partial modes are made costlier        
    otherwise
        warning('Unexpected slip type.\n {Default perfect slip is applied');        
 end
        
% criteria = 'minplasticspin';

gmatrix = gmatrix0;
lg =  length(gmatrix);
e_ext= [1,0,0;0,0,0;0,0,-1];
b = zeros(5,lg);
n_steps = max_strain/incr_e;
new_gmatrix(:,1) = gmatrix;

%% Calculation Block

for j=1:1:n_steps

for i=1:1:lg 
    DCMtrx = matrix(gmatrix(i));  
    e_grain=DCMtrx'*e_ext*DCMtrx; 
    b = [e_grain(1,1);e_grain(2,2);2*e_grain(2,3);2*e_grain(1,3);2*e_grain(1,2)]; 
    AllSolutions = calcLPSolution(b,c); 
    [spin, UniqueSolution] = calculate_spin(AllSolutions,criteria);
    spin = incr_e*spin;
    gmatrix(i) = applySpin(-spin,DCMtrx); 
    new_gmatrix(i,j+1) = gmatrix(i);
%     ActiveSS(i,j).ss= UniqueSolution.B;
%     ActiveSS(i,j).gamma = UniqueSolution.xb;
end
end
%%

p = Euler(new_gmatrix(:,(j+1)))/degree;
delp = p - p0;


%%

toc