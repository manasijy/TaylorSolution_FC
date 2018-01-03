function TaylorFCDefTexSim()

tic
load('FCC_SS24_Set.mat','SlipSystem');

%%
prompt = 'The euler angle file name with .txt extension \n';
g_file = input(prompt);                      
g = fopen(g_file);   
g_matrix = textscan(g, '%f %f %f'); 
fclose(g);

%% Storing all orientation as Euler angle matrix and also as direction cosines

gmatrixR =[g_matrix{1,1},g_matrix{1,2},g_matrix{1,3}]*degree;% changed to radians
dcm_gmatrix = angle2dcm(gmatrixR(:,1),gmatrixR(:,2),gmatrixR(:,3),'ZXZ');% it is a 3x3xlg size matrix
lg =  length(gmatrixR);

%% Initialization of variables

e_ext= [1,0,0;0,0,0;0,0,-1];
b = zeros(5,1);%lg);
max_strain = 1;
incr_e = 0.1; % incremental strain per step
n_steps = max_strain/incr_e;
% newdcm_gmatrix = dcm_gmatrix;
new_gmatrixR = zeros(lg,n_steps+1,3); %The output matrix with updated orientation after each deformation step
for j=1:1:n_steps+1
    if j==1, new_gmatrixR(:,1,:) = gmatrixR; continue; end
for i=1:1:lg, 
    e_grain=dcm_gmatrix(:,:,i)'*e_ext*dcm_gmatrix(:,:,i);
    b = [e_grain(1,1);e_grain(2,2);2*e_grain(2,3);2*e_grain(1,3);2*e_grain(1,2)]; 
    spin = calcPlasticSpin(b);
    spin_min = incr_e*spin.min;
    dcm_gmatrix(:,:,i) = rotategrain(spin_min)*dcm_gmatrix(:,:,i);
    [phi1, phi, phi2] = dcm2angle(dcm_gmatrix(:,:,i),'ZXZ');
    new_gmatrixR(i,j,:) = [phi1, phi, phi2];
end
end
% new_gmatrixR = new_gmatrixR/degree;
save('OutputTexture.mat', 'new_gmatrixR')
toc
end


function [rotmatrix] = rotategrain(spin)

rotangle = sqrt(spin(1,2)^2 + spin(1,3)^2 + spin(2,3)^2);
rotaxis = [-spin(2,3)/rotangle, -spin(3,1)/rotangle, -spin(1,2)/rotangle];
r = [rotaxis rotangle];
rotmatrix = vrrotvec2mat(r);
end

function [spin] = calcPlasticSpin(b)
load('FCC_SS24_Set.mat','SlipSystem');
A = [0,-0.408248290463863,-0.408248290463863,-0.408248290463863,-0.408248290463863,0,0.408248290463863,...
    -0.408248290463863,0,-0.408248290463863,0.408248290463863,0,...
    0,0.408248290463863,0.408248290463863,0.408248290463863,0.408248290463863,0,-0.408248290463863,...
    0.408248290463863,0,0.408248290463863,-0.408248290463863,0;
    
 -0.408248290463863,0,0.408248290463863,0.408248290463863,0,...
 -0.408248290463863,-0.408248290463863,0,-0.408248290463863,0.408248290463863,...
 0,0.408248290463863,0.408248290463863,0,-0.408248290463863,...
 -0.408248290463863,0,0.408248290463863,0.408248290463863,0,...
 0.408248290463863,-0.408248290463863,0,-0.408248290463863;
 
 0,0.408248290463863,0.408248290463863,0.408248290463863,0.408248290463863,0,...
 0.408248290463863,-0.408248290463863,0,-0.408248290463863,0.408248290463863,...
 0,0,-0.408248290463863,-0.408248290463863,-0.408248290463863,-0.408248290463863,...
 0,-0.408248290463863,0.408248290463863,0,0.408248290463863,...
 -0.408248290463863,0;
 
 0.408248290463863,0,-0.408248290463863,0.408248290463863,0,...
 -0.408248290463863,0.408248290463863,0,0.408248290463863,0.408248290463863,...
 0,0.408248290463863,-0.408248290463863,0,0.408248290463863,...
 -0.408248290463863,0,0.408248290463863,-0.408248290463863,0,...
 -0.408248290463863,-0.408248290463863,0,-0.408248290463863;
 
 -0.408248290463863,-0.408248290463863,0,0,0.408248290463863,...
 0.408248290463863,0,0.408248290463863,0.408248290463863,0,0.408248290463863,0.408248290463863,...
 0.408248290463863,0.408248290463863,0,0,-0.408248290463863,...
 -0.408248290463863,0,-0.408248290463863,-0.408248290463863,0,-0.408248290463863,-0.408248290463863];
c = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
k = 0;
tol=1e-10; 
maxit=50; 

%%
[m,n]=size(A); b=b(:); c=c(:); it=0; 
if (length(c)~=n || length(b)~=m),error('wrong dimensions'); end
D=sign(sign(b)+.5); 
D = diag(D);                    % initial (inverse) basis matrix
A = [A D];                      % incorporate slack/artificial variables.Utilizing 2 phase method

%%
B = n+1:n+m;                    
N = 1:n;                        
phase=1; xb=abs(b); s=[zeros(n,1);ones(m,1)];   % supercost 

%%
while phase<3,
   df=-1; t=inf;
   yb= D'*s(B);  
   while (it < maxit)
      if isempty(N), break, end     % no freedom for minimization
      r = s(N) - [A(:,N)]'*yb;      % Transpose both sides to get the standard form of this equation
      [rmin,q] = min(r);            % determine minimum of reduced cost % if reduced cost is zero,
                                    %there is no improvement i.e. another optimal solution
      if rmin>=-tol*(norm(s(N),inf)+1), break, end % optimal!            
      it=it+1;
      if df>=0                      % apply Bland's rule to avoid cycling
         if maxit==inf,
            disp(['LINPROG(',int2str(it),'): warning! degenerate vertex']);
         end
         J=find(r<0); Nq=min(N(J)); q=find(N==Nq); % This line gets the first most column in NB which gives a -ve reduced cost. 
                                    % q is the location of that column in N                                                                                                     
      end 
      d = D*A(:,N(q)); % Descent Direction: d = -inv(B)*Nq.If dq >= 0, then the linear program is unbounded. Here d<=0 =>unbounded
      I=find(d>tol); %find out the entry no of those d values which are -ve
      if isempty(I), disp('Solution is unbounded'); it=-it; break; end
      xbd=xb(I)./d(I); [r,p]=min(xbd); p=I(p); % r is the step length alpha, p is position of that alpha, p will be going out
      if df>=0,                     % apply Bland's rule to avoid cycling
         J=find(xbd==r); Bp=min(B(I(J))); p=find(B==Bp); 
      end 
      xb= xb - r*d; xb(p)=r;        % update x, bsic var xb(p) is replaced with alpha 
      df=r*rmin;                    % change in f 
      v = D(p,:)/d(p);              % row vector
      yb= yb + v'*( s(N(q)) - d'*s(B) );
      d(p)=d(p)-1;
      D = D - d*v;                  % update inverse basis matrix
      t=B(p); B(p)=N(q); % pth column of B is replaced with qth column of N. q is the one with minimum reduced cost. it has come in. 
      if t>n+k, N(q)=[]; else N(q)=t; end %The pth column of B is brought into N.
                                          %If t was artificial var then one col of N is reduced:single phase
   end                             % end of phase
   xb=xb+D*(b-A(:,B)*xb);          % iterative refinement
   I=find(xb<0);                   % must be due to rounding error
   if I, xb(I)=xb(I)-xb(I); end    % so correct
   if phase==2 || it<0, break; end; % B, xb,n,m,res=A(:,B)*xb-b %| changed to ||
   if xb'*s(B)>tol,it=-it; disp('no feasible solution'); break;  end % if 
   phase=phase+1;                  % Go to Phase 2
   s=1e6*norm(c,'inf')*s; s(1:n)=c;% tol=tol*norm(s,inf);
end
x=sparse(n,1); x(B)=xb; x=x(1:n); if n<21, x=full(x); end
fval=c'*x;
if it>=maxit, disp('too many iterations'); it=-it; end

OP.B = B; OP.N = N; OP.D = D; OP.r = r; OP.xb = xb; OP.yb = yb; OP.var = var(xb); 
        
%% Getting multiple solutions and passing on the plastic spin as output

TaylorSolution = treeSol_function(OP,b);
n_ts = numel(TaylorSolution);
varArray =  cell2mat({TaylorSolution.var});
minVar = find(varArray == min(varArray)); 
maxVar = find(varArray == max(varArray));
minVarSolution = TaylorSolution(minVar(1)); % Need to check if two results are also possible
maxVarSolution = TaylorSolution(maxVar(1));
% n_mvs = numel(minVarSolution.B);
% n_mxvs = numel(maxVarSolution.B);I think both will be 5
SSmin= SlipSystem(minVarSolution.B);
SSmax= SlipSystem(maxVarSolution.B);
shearmin = minVarSolution.xb;
shearmax = maxVarSolution.xb;
spin.min = zeros(3,3);
spin.max = zeros(3,3);

for ii=1:1:5, 
    spin.min = spin.min + shearmin(ii)*SSmin(ii).q.M;
    spin.max = spin.max + shearmax(ii)*SSmax(ii).q.M;
end
end

%%

function [sol] = treeSol_function(stPt,b)

new(1)=stPt;
sol(1)=stPt;

while ~isempty(new)
    first = new(1);
    new(1) = [];
    temp = multipleSol_function(first,b);

    for i=1:1:numel(temp)
        repeat = 0;
       
            for j = 1:1:numel(sol)
            if isequaln(sort(sol(j).B),sort(temp(i).B))
                repeat=1;
                break;
            end
            end
        
        if ~repeat                 
           new = [new,temp(i)];
           sol = [sol,temp(i)];
        end
    end
end
end

%%

function [slip] = multipleSol_function(stPt,b)

%%
load('FCC_SS24_Set.mat','SlipSystem');
A = [0,-0.408248290463863,-0.408248290463863,-0.408248290463863,-0.408248290463863,0,0.408248290463863,...
    -0.408248290463863,0,-0.408248290463863,0.408248290463863,0,...
    0,0.408248290463863,0.408248290463863,0.408248290463863,0.408248290463863,0,-0.408248290463863,...
    0.408248290463863,0,0.408248290463863,-0.408248290463863,0;
    
 -0.408248290463863,0,0.408248290463863,0.408248290463863,0,...
 -0.408248290463863,-0.408248290463863,0,-0.408248290463863,0.408248290463863,...
 0,0.408248290463863,0.408248290463863,0,-0.408248290463863,...
 -0.408248290463863,0,0.408248290463863,0.408248290463863,0,...
 0.408248290463863,-0.408248290463863,0,-0.408248290463863;
 
 0,0.408248290463863,0.408248290463863,0.408248290463863,0.408248290463863,0,...
 0.408248290463863,-0.408248290463863,0,-0.408248290463863,0.408248290463863,...
 0,0,-0.408248290463863,-0.408248290463863,-0.408248290463863,-0.408248290463863,...
 0,-0.408248290463863,0.408248290463863,0,0.408248290463863,...
 -0.408248290463863,0;
 
 0.408248290463863,0,-0.408248290463863,0.408248290463863,0,...
 -0.408248290463863,0.408248290463863,0,0.408248290463863,0.408248290463863,...
 0,0.408248290463863,-0.408248290463863,0,0.408248290463863,...
 -0.408248290463863,0,0.408248290463863,-0.408248290463863,0,...
 -0.408248290463863,-0.408248290463863,0,-0.408248290463863;
 
 -0.408248290463863,-0.408248290463863,0,0,0.408248290463863,...
 0.408248290463863,0,0.408248290463863,0.408248290463863,0,0.408248290463863,0.408248290463863,...
 0.408248290463863,0.408248290463863,0,0,-0.408248290463863,...
 -0.408248290463863,0,-0.408248290463863,-0.408248290463863,0,-0.408248290463863,-0.408248290463863];

r = stPt.r;
B = stPt.B;
N = stPt.N;
D = stPt.D;
xb = stPt.xb;
yb = stPt.yb;
s = [1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1000000;1000000;1000000;1000000;1000000];
slip(1) = stPt;

%%
zero_r = find(abs(r)<1e-10);
if ~isempty(zero_r)
jj = 1;
l_s = numel(zero_r);                                             
zero_d = D*A(:,N(zero_r)); % = inv(B)*N(q)
round_d = round(zero_d);
all1 =  ones(5,1);
xzero = find(xb<0.0001);
all1(xzero) = 0;

for i = 1:1:l_s
        B1 = B;
        N1 = N;
        mark = all1-round_d(:,i); %Check to remove cases where x becomes -ve
        if (min(mark) < 0), continue; end
        
        rp = find(round_d(:,i) == 1);
        [xmin,xminp] = min(xb(rp));       
         p1 = rp(xminp); q1 = zero_r(i);d1 = zero_d(:,i);
         % p1 is leaving and q1 is entering
        v1 = D(p1,:)/d1(p1);             
      yb1= yb + v1'*( s(N1(q1)) - d1'*s(B1) );
      d1(p1)=d1(p1)-1;
      D1 = D - d1*v1;                          
        t1 = B1(p1);
        B1(p1) = N1(q1);
        N1(q1) = t1;        
        r1 = s(N1) - [A(:,N1)]'*yb1;                   
        xb1= xb - xmin*d1; xb1(p1)=xmin;
        xb1=xb1+D1*(b-A(:,B1)*xb1);          % iterative refinement
        I=find(xb1<0);                   % must be due to rounding error
        if I, xb1(I)=xb1(I)-xb1(I); end    % so correct
              
        slip(jj).B = B1; slip(jj).N = N1; slip(jj).D = D1; slip(jj).r = r1;
        slip(jj).xb = xb1; slip(jj).yb = yb1; slip(jj).var = var(xb1);
        
        jj = jj+1;
end
end
end