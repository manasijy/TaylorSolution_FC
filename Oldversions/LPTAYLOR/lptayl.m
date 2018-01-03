
clear

%% Input
fprintf('This program calculates texture evolution with deformation strain \nfor e11=1 and e33=-1 rolling \n');
prompt = 'The euler angle file name without .txt extension \n';
g_file = [input(prompt), '.txt']; 
prompt = 'Please enter maximum strain e.g. 2,3 etc \n';
max_strain = input(prompt);
prompt = 'Please enter incremental strain e.g. 0.05,0.1 etc \n';
incr_e = input(prompt);
tic
g = fopen(g_file);   
g_matrix = textscan(g, '%f %f %f'); 
fclose(g);

%% Initialization
gmatrix0 = rotation('Euler',[g_matrix{1,1},g_matrix{1,2},g_matrix{1,3}]*degree);
gmatrix = gmatrix0;
lg =  length(gmatrix);
e_ext= [1,0,0;0,0,0;0,0,-1];
b = zeros(5,lg);
% max_strain = 3;
% incr_e = 0.1;                   % incremental strain per step
n_steps = max_strain/incr_e;
new_gmatrix = gmatrix;

%%
for j=1:1:n_steps+1
    if j==1, new_gmatrix(:,j) = gmatrix; continue; end
for i=1:1:lg, 
    DCMtrx = (matrix(gmatrix(i)))';  
    e_grain=DCMtrx'*e_ext*DCMtrx;
    b = [e_grain(1,1);e_grain(2,2);2*e_grain(2,3);2*e_grain(1,3);2*e_grain(1,2)]; 
    [spin, solution] = calcPlasticSpin(b); 
    spin_min = incr_e*spin.min; 
    gmatrix(i) = applySpin(spin_min,DCMtrx);    
    new_gmatrix(i,j) = gmatrix(i);
    
    ActiveSS(i,j-1).ss= solution.B;
    ActiveSS(i,j-1).gamma = solution.xb;
end
end

%% Output

folderPath = resmet_input(g_file,gmatrix);
eu_gmat = Euler(new_gmatrix)/degree;
eu_filename = [folderPath '\' 'eulerAngles.txt'];
save(eu_filename,eu_gmat,'-ascii');
ActiveSSData = [folderPath '\' 'ActiveSSData.mat'];
save(ActiveSSData,'ActiveSS');
toc

%  s =0.05*sspin.M
%  rotdcmatrix = spin2mat(s)
%  dcg1=  rotdcmatrix'*dcg0
%  [r1 r2 r3] = dcm2angle(dcg1', 'ZXZ')
