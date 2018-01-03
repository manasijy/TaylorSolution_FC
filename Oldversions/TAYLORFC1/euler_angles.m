function  euler_angles(name,new_gmatrix)

%%
%This fuction takes the output of main_LP.. function and creates a folder
%which contains resmet compatible input file consisting of euler angles at
%each deformation step
%%

[r, c] = size(new_gmatrix);
foldername = [name 'EU' date];
if exist(foldername,'dir'), foldername = [foldername  num2str(random('Poisson',1))]; end
    
mkdir(foldername);
path = fullfile(cd,foldername); %cd reads the current path in string format
for j=1:1:c
    eu_gmat = Euler(new_gmatrix(:,j))/degree;
    eu_filename = [path '\' num2str(j) 'eulerAngles.txt'];
    save(eu_filename,'eu_gmat','-ascii');

end
end