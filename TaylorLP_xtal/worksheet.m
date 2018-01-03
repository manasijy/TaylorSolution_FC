%  testd =(TaylorSolution(2).D)*A(:,TaylorSolution(2).N([8,14,15]));
clear;
load('SlipSystem24.mat','SlipSystem'); 
A = [[SlipSystem(1).m.M(1,1) , SlipSystem(2).m.M(1,1) , SlipSystem(3).m.M(1,1)...
     , SlipSystem(4).m.M(1,1),SlipSystem(5).m.M(1,1),SlipSystem(6).m.M(1,1)...
     ,SlipSystem(7).m.M(1,1),SlipSystem(8).m.M(1,1),SlipSystem(9).m.M(1,1)...
     ,SlipSystem(10).m.M(1,1),SlipSystem(11).m.M(1,1),SlipSystem(12).m.M(1,1)...
     ,SlipSystem(13).m.M(1,1) , SlipSystem(14).m.M(1,1) , SlipSystem(15).m.M(1,1)...
     , SlipSystem(16).m.M(1,1),SlipSystem(17).m.M(1,1),SlipSystem(18).m.M(1,1)...
     ,SlipSystem(19).m.M(1,1),SlipSystem(20).m.M(1,1),SlipSystem(21).m.M(1,1)...
     ,SlipSystem(22).m.M(1,1),SlipSystem(23).m.M(1,1),SlipSystem(24).m.M(1,1);
 SlipSystem(1).m.M(2,2) , SlipSystem(2).m.M(2,2) , SlipSystem(3).m.M(2,2)...
     , SlipSystem(4).m.M(2,2),SlipSystem(5).m.M(2,2),SlipSystem(6).m.M(2,2)...
     ,SlipSystem(7).m.M(2,2),SlipSystem(8).m.M(2,2),SlipSystem(9).m.M(2,2)...
     ,SlipSystem(10).m.M(2,2),SlipSystem(11).m.M(2,2),SlipSystem(12).m.M(2,2)...
     SlipSystem(13).m.M(2,2) , SlipSystem(14).m.M(2,2) , SlipSystem(15).m.M(2,2)...
     , SlipSystem(16).m.M(2,2),SlipSystem(17).m.M(2,2),SlipSystem(18).m.M(2,2)...
     ,SlipSystem(19).m.M(2,2),SlipSystem(20).m.M(2,2),SlipSystem(21).m.M(2,2)...
     ,SlipSystem(22).m.M(2,2),SlipSystem(23).m.M(2,2),SlipSystem(24).m.M(2,2)];
 2*[SlipSystem(1).m.M(2,3) , SlipSystem(2).m.M(2,3) , SlipSystem(3).m.M(2,3)...
     , SlipSystem(4).m.M(2,3),SlipSystem(5).m.M(2,3),SlipSystem(6).m.M(2,3)...
     ,SlipSystem(7).m.M(2,3),SlipSystem(8).m.M(2,3),SlipSystem(9).m.M(2,3)...
     ,SlipSystem(10).m.M(2,3),SlipSystem(11).m.M(2,3),SlipSystem(12).m.M(2,3)...
     ,SlipSystem(13).m.M(2,3) , SlipSystem(14).m.M(2,3) , SlipSystem(15).m.M(2,3)...
     , SlipSystem(16).m.M(2,3),SlipSystem(17).m.M(2,3),SlipSystem(18).m.M(2,3)...
     ,SlipSystem(19).m.M(2,3),SlipSystem(20).m.M(2,3),SlipSystem(21).m.M(2,3)...
     ,SlipSystem(22).m.M(2,3),SlipSystem(23).m.M(2,3),SlipSystem(24).m.M(2,3);
SlipSystem(1).m.M(1,3) , SlipSystem(2).m.M(1,3) , SlipSystem(3).m.M(1,3)...
     , SlipSystem(4).m.M(1,3),SlipSystem(5).m.M(1,3),SlipSystem(6).m.M(1,3)...
     ,SlipSystem(7).m.M(1,3),SlipSystem(8).m.M(1,3),SlipSystem(9).m.M(1,3)...
     ,SlipSystem(10).m.M(1,3),SlipSystem(11).m.M(1,3),SlipSystem(12).m.M(1,3)...
     ,SlipSystem(13).m.M(1,3) , SlipSystem(14).m.M(1,3) , SlipSystem(15).m.M(1,3)...
     , SlipSystem(16).m.M(1,3),SlipSystem(17).m.M(1,3),SlipSystem(18).m.M(1,3)...
     ,SlipSystem(19).m.M(1,3),SlipSystem(20).m.M(1,3),SlipSystem(21).m.M(1,3)...
     ,SlipSystem(22).m.M(1,3),SlipSystem(23).m.M(1,3),SlipSystem(24).m.M(1,3);
SlipSystem(1).m.M(1,2) , SlipSystem(2).m.M(1,2) , SlipSystem(3).m.M(1,2)...
     , SlipSystem(4).m.M(1,2),SlipSystem(5).m.M(1,2),SlipSystem(6).m.M(1,2)...
     ,SlipSystem(7).m.M(1,2),SlipSystem(8).m.M(1,2),SlipSystem(9).m.M(1,2)...
     ,SlipSystem(10).m.M(1,2),SlipSystem(11).m.M(1,2),SlipSystem(12).m.M(1,2)...
     ,SlipSystem(13).m.M(1,2) , SlipSystem(14).m.M(1,2) , SlipSystem(15).m.M(1,2)...
     , SlipSystem(16).m.M(1,2),SlipSystem(17).m.M(1,2),SlipSystem(18).m.M(1,2)...
     ,SlipSystem(19).m.M(1,2),SlipSystem(20).m.M(1,2),SlipSystem(21).m.M(1,2)...
     ,SlipSystem(22).m.M(1,2),SlipSystem(23).m.M(1,2),SlipSystem(24).m.M(1,2);... 
    ]];


r = [2.00000000000000;2.00000000000000;2;2.00000000000000;2.00000000000000;...
    1.00000000000000;2.00000000000000;-4.44089209850063e-16;1;1.00000000000000;...
    1.00000000000000;0.999999999999999;1.00000000000000;0;4.44089209850063e-16;1.00000000000000;2;1;2];
B = [19,12,17,18,13];
N = [1,2,5,6,7,8,9,14,11,4,3,15,16,10,21,20,22,23,24];
D = [-1.22474487139159,1.66533453693773e-16,-1.22474487139159,8.27273252079442e-18,-1.58260721172979e-16;...
    2.22044604925031e-16,-5.55111512312578e-17,-1.22474487139159,1.22474487139159,1.22474487139159;...
    1.22474487139159,1.91485146782310e-16,-1.22474487139159,6.66133814775094e-16,0;...
    0,1.22474487139159,1.22474487139159,-4.44089209850063e-16,-1.22474487139159;...
    1.22474487139159,1.22474487139159,1.22474487139159,-1.22474487139159,4.62100654151109e-16];
xb = [0.0186064918682347;0.0161110278532987;0.0184651536349038;0.985040311823264;1.41112528465517];
yb = [1.22474487139159;2.44948974278318;-1.22474487139159;6.29966591578187e-16;4.07921986653155e-16];
s = [1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1000000;1000000;1000000;1000000;1000000];
b = [-5.77010921345182e-05;0.992403876506104;-0.0151344359013386;-0.174967230668978;0.172987393925089];


%%
zero_r = find(abs(r)<1e-10);
if ~isempty(zero_r)
jj = 0;
l_s = numel(zero_r);                                             
zero_d = D*A(:,N(zero_r)); % = inv(B)*N(q)
round_d = round(zero_d);
all1 =  ones(5,1);
xzero = find(xb<0.0001);
all1(xzero) = 0;

for i = 1:1:l_s
%         i=3;
jj=jj+1;
        mark = all1-round_d(:,i); %Check to remove cases where x becomes -ve
        if (min(mark) < 0), continue; end
        
        rp = find(round_d(:,i) == 1);
        [xmin,xminp] = min(xb(rp));       
         p1 = rp(xminp); q1 = zero_r(i);d1 = zero_d(:,i);
         % p1 is leaving and q1 is entering
        v1 = D(p1,:)/d1(p1);             
      yb1= yb + v1'*( s(N(q1)) - d1'*s(B) );
      d1(p1)=d1(p1)-1;
      D1 = D - d1*v1;                          
        t1 = B(p1);
        B(p1) = N(q1);
        N(q1) = t1;        
        r1 = s(N) - [A(:,N)]'*yb1;                   
        xb1= xb - xmin*d1; xb1(p1)=xmin;
        xb1=xb1+D1*(b-A(:,B)*xb1);          % iterative refinement
        I=find(xb1<0);                   % must be due to rounding error
        if I, xb1(I)=xb1(I)-xb1(I); end    % so correct
              
        slip(jj).B = B; slip(jj).N = N; slip(jj).D = D1; slip(jj).r = r1;
        slip(jj).xb = xb1; slip(jj).yb = yb1; slip(jj).M = sum(xb1);slip(jj).s=s;
        
%         jj = jj+1;
end
end
% end