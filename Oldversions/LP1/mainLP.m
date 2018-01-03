tic
clear

%%
% gvector = [20,45,50];
% e_ext =[1,0,0;0,0,0;0,0,-1];

prompt = 'The euler angle file name with .txt extension \n';
g_vectorfile = input(prompt);                      
g = fopen(g_vectorfile);   
g_matrix = textscan(g, '%f %f %f'); 
fclose(g);
lg =  length(g_matrix{1,1});

%%
e_ext= [1,0,0;0,-1,0;0,0,0];


        
%% for one g vector


OP = LP_function(g,e_ext);
TaylorSolution = treeSol_function(OP,A,b);

toc