
clear

%% Input Block
fprintf('This program calculates texture evolution with deformation strain \nfor e11=1 and e33=-1 rolling \n');
prompt = 'The euler angle file name without .txt extension \n';
name = input(prompt);
g_file = [name, '.txt']; 
prompt = 'Please enter maximum strain e.g. 2,3 etc \n';
max_strain = input(prompt);
prompt = 'Please enter incremental strain e.g. 0.05,0.1 etc \n';
incr_e = input(prompt);
tic
g = fopen(g_file);   
g_matrix = textscan(g, '%f %f %f'); 
fclose(g);

%% Initialization

% Criteria will be given as input in later development. But for present
% minvar is taken as the single unique solution selection criteria. Others
% are 'maxvar' and 'minplasticspin'
criteria = 'minvar';
% criteria = 'maxvar';
% criteria = 'minplasticspin';

gmatrix0 = rotation('Euler',[g_matrix{1,1},g_matrix{1,2},g_matrix{1,3}]*degree);
gmatrix = gmatrix0;
lg =  length(gmatrix);
e_ext= [1,0,0;0,0,0;0,0,-1];
b = zeros(5,lg);
n_steps = max_strain/incr_e;
new_gmatrix(:,1) = gmatrix;
%% Calculation Block

for j=1:1:n_steps

for i=1:1:lg, 
    DCMtrx = matrix(gmatrix(i));  
    e_grain=DCMtrx'*e_ext*DCMtrx; 
    b = [e_grain(1,1);e_grain(2,2);2*e_grain(2,3);2*e_grain(1,3);2*e_grain(1,2)]; 
    AllSolutions = calcLPSolution(b); 
    [spin, UniqueSolution] = calculate_spin(AllSolutions,criteria);
    spin = incr_e*spin;
    gmatrix(i) = applySpin(-spin,DCMtrx); 
    new_gmatrix(i,j+1) = gmatrix(i);
    ActiveSS(i,j).ss= UniqueSolution.B;
    ActiveSS(i,j).gamma = UniqueSolution.xb;
end
end

%% Output Block
prompt = 'Do you want to save results; yes or no \n';
reply = input(prompt,'s');
if strcmp(reply,'yes')
    folderpath = resmet_input(name,new_gmatrix);

    ActiveSSData = [folderpath '\' 'ActiveSSData.mat'];
    save(ActiveSSData,'ActiveSS');
    
    %%%Following lines to save o/p euler angles in file
    %     eu_gmat = Euler(new_gmatrix)/degree;
    %     euler_angles(name,new_gmatrix);
else break
end
%%

toc

