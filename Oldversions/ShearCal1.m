%% This program creates slip system set including schimd factor for each system
% Both positive and negative SS are taken in to consideration


e11 = 1; e22 = -1; e23 = 0; e13 = 0; e12 = 0;
e_ext = [e11;e22;2*e23;2*e13;2*e12];
B = e_ext;
sol = struct('no',zeros(1,5),'shear',zeros(1,5));
CS = crystalSymmetry('-43m');

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
SlipSystem(10).n=[1,1,-1];SlipSystem(10).b =[-1,1,0];
SlipSystem(11).n=[1,1,-1]; SlipSystem(11).b=[1,0,1];
SlipSystem(12).n=[1,1,-1]; SlipSystem(12).b =[0,1,1];


for i = 1:1:12;

n = SlipSystem(i).n;
b = SlipSystem(i).b;
hkl = Miller(n(1),n(2),n(3),CS);
uvw = Miller(b(1),b(2),b(3),CS);

R = SchmidTensor(hkl,uvw);  %,'generalized');
% Rmatrix = matrix(R);
SlipSystem(i).SchmidT = R; %.M;

end

sol_set = 1;
count = 0;
for i = 1:1:5
    for j = (i+1):1:8

        for k = (j+1):1:9

            for l = (k+1):1:11

                for m = (l+1):1:12
                        

                        A(1,:) = [SlipSystem(1,i).SchmidT.M(1,1) , SlipSystem(1,j).SchmidT.M(1,1) , SlipSystem(1,k).SchmidT.M(1,1)...
                            , SlipSystem(1,l).SchmidT.M(1,1),SlipSystem(1,m).SchmidT.M(1,1)];
                        A(2,:) = [SlipSystem(1,i).SchmidT.M(2,2) , SlipSystem(1,j).SchmidT.M(2,2) , SlipSystem(1,k).SchmidT.M(2,2)...
                            , SlipSystem(1,l).SchmidT.M(2,2),SlipSystem(1,m).SchmidT.M(2,2)];
                        A(3,:) = 2*[SlipSystem(1,i).SchmidT.M(2,3) , SlipSystem(1,j).SchmidT.M(2,3) , SlipSystem(1,k).SchmidT.M(2,3)...
                            , SlipSystem(1,l).SchmidT.M(2,3),SlipSystem(1,m).SchmidT.M(2,3)];
                        A(4,:) = 2*[SlipSystem(1,i).SchmidT.M(1,3) , SlipSystem(1,j).SchmidT.M(1,3) , SlipSystem(1,k).SchmidT.M(1,3)...
                            , SlipSystem(1,l).SchmidT.M(1,3),SlipSystem(1,m).SchmidT.M(1,3)];
                        A(5,:) = 2*[SlipSystem(1,i).SchmidT.M(1,2) , SlipSystem(1,j).SchmidT.M(1,2) , SlipSystem(1,k).SchmidT.M(1,2)...
                            , SlipSystem(1,l).SchmidT.M(1,2),SlipSystem(1,m).SchmidT.M(1,2)];
                        if (det(A)>0.001)
                        X_G = A\B; %for AX=B type eq  
                                              
                        sol(sol_set).no = [i,j,k,l,m];
                        sol(sol_set).shear = X_G';
                        sol(sol_set).total = sum(abs(X_G));
                        sol_set = sol_set+1;
                        end
                        count = count+1;
                end
           end
         end
    end
end

save('shearG2.mat','sol');
fid= fopen('shearfile1.txt','a+');
fprintf(fid,' Slip system no \n Shear in corresponding slip system\n\n');

for n_s=1:1:sol_set-1
    nx = (sol(n_s).no);
    ns = (sol(n_s).shear);
    fprintf(fid,'Slip System Set %d \n',n_s);
    fprintf(fid,'%d %d %d  %d  %d\n',nx);
    fprintf(fid,'%+3.2f %+3.2f %+3.2f %+3.2f %+3.2f Total Shear = %3.2f \n\n',ns,sol(n_s).total);
end
