
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
lg =  length(g_matrix);
e_ext= [1,0,0;0,0,0;0,0,-1];
b = zeros(5,lg);
n_steps = max_strain/incr_e;


%% Calculation Block

gmatrix = euler2gmatrix([g_matrix{1,1},g_matrix{1,2},g_matrix{1,3}]*degree);

for j=1:1:n_steps+1
    if j==1, new_gmatrix = gmatrix; continue; end
    
for i=1:1:lg,  
    DCMtrx = DC(new_gmatrix(i));
    e_grain=DCMtrx'*e_ext*DCMtrx; % the component tr matrix is inverse(transpose) of axis tr dc matrix
    b = [e_grain(1,1);e_grain(2,2);2*e_grain(2,3);2*e_grain(1,3);2*e_grain(1,2)]; 
    AllSolutions = calcLPSolution(b); 
    [spin, UniqueSolution] = calculate_spin(AllSolutions,criteria); %%%%%%!!!!!
    spin = incr_e*spin;
    gmatrix(i) = applySpin(spin,DCMtrx);   %applySpin(-spin,DCMtrx');  
    new_gmatrix(i,j) = gmatrix(i);
    ActiveSS(i,j-1).ss= UniqueSolution.B;
    ActiveSS(i,j-1).gamma = UniqueSolution.xb;
end
end

% %% Output Block
% prompt = 'Do you want to save results; yes or no \n';
% reply = input(prompt,'s');
% if strcmp(reply,'yes')
%     folderpath = resmet_input(name,new_gmatrix);
%     eu_gmat = Euler(new_gmatrix)/degree;
%     ActiveSSData = [folderpath '\' 'ActiveSSData.mat'];
%     save(ActiveSSData,'ActiveSS');
%     euler_angles(name,new_gmatrix);
% else break
% end
%%

toc

