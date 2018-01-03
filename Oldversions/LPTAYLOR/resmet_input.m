function [path] = resmet_input(g_file,gmatrix)

%%
%This fuction takes the output of main_LP.. function and creates a folder
%which contains resmet compatible input file consisting of euler angles at
%each deformation step
%%

[r, c] = size(gmatrix);
foldername = [g_file date];
% if exist(foldername,'dir')
%     prompt 'nputdate already exists input some indentification for
mkdir(foldername);
% path = [cd '\ResmetInput' date]; %cd reads the current path in string format
path = fullfile(cd,foldername); %cd reads the current path in string format
for j=1:1:c
    filename = ['Def_step' num2str(j-1) '.txt'];
    fpath = fullfile(path,filename);
    fid= fopen(fpath,'w');
    fprintf(fid,'%s \n',num2str(r));
    fclose(fid);
        for i=1:1:r
            fid= fopen(fpath,'a+');
            fprintf(fid,'\t%4.2f\t%4.2f\t%4.2f\t%d\n',Euler(gmatrix(i,j))/degree,1);
            if i==r, fclose(fid); end
        end
    
end
    
    