tic
clear

%%
gvector = [90,90,90];
e_ext =[1,0,0;0,0,0;0,0,-1];

%%

c = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
k = 0;
tol=1e-10; 
maxit=50; 
DCMtrx = DC_matrix_function(gvector(1),gvector(2),gvector(3));
[e_grain]= transform_e_function(e_ext,DCMtrx);
b = [e_grain(1,1);e_grain(2,2);2*e_grain(2,3);2*e_grain(1,3);2*e_grain(1,2)]; 

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

OP.B = B; OP.N = N; OP.D = D; OP.r = r;
        OP.xb = xb; OP.yb = yb; OP.M = sum(xb);OP.s=s;
        
%%

TaylorSolution = treeSol_function(OP,A,b);

%%
n_ts = numel(TaylorSolution);
for i=1:n_ts, TaylorSolution(i).var = var(TaylorSolution(i).xb); end
varArray =  cell2mat({TaylorSolution.var});
minVar = find(varArray == min(varArray));
maxVar = find(varArray == max(varArray));
minVarSolution = TaylorSolution(minVar(1)); % Need to check if two results are also possible
maxVarSolution = TaylorSolution(maxVar(1));
n_mvs = numel(minVarSolution.B);
n_mvs1 = numel(maxVarSolution.B);
SSmin= SlipSystem(minVarSolution.B);
SSmax= SlipSystem(maxVarSolution.B);
shearmin = minVarSolution.xb;
shearmax = maxVarSolution.xb;
Minspin = zeros(3,3);
Maxspin = Minspin;
cs = crystalSymmetry('-43m');
ss = crystalSymmetry('mmm');

%%
o0 = orientation('Euler',gvector*degree,cs,ss);
plotPDF(o0,Miller(1,1,1,cs),'complete');
%% Min spin
for ii=1:n_mvs, 
    Minspin = Minspin + shearmin(ii)*SSmin(ii).q.M;
end
% spin = Minspin;
rotmin = rotation('matrix',-Minspin);
o1=orientation(rotmin,cs,ss)*o0;
figure
plotPDF(o1,Miller(1,1,1,cs),'complete');

%% Max spin

for ii=1:n_mvs1, 
    Maxspin = Maxspin + shearmax(ii)*SSmax(ii).q.M;
end
% spin = Minspin;
rotmax = rotation('matrix',-Maxspin);
o2=orientation(rotmax,cs,ss)*o0;
figure
plotPDF(o2,Miller(1,1,1,cs),'complete');

Euler(o0)
Euler(o1)
Euler(o2)

%
cs = crystalSymmetry('-43m');
ss = crystalSymmetry('mmm');
name = {'brass','goss','copper'};
for i=1:3
    figure
    o = orientation(name{i},cs,ss);
    plotPDF(o,Miller(1,1,2,cs),'complete','FigureTitle',name{i});
end
%%

function [sol] = treeSol_function(stPt,A,b)

new(1)=stPt;
sol(1)=stPt;

while ~isempty(new)
    first = new(1);
    new(1) = [];
    temp = multipleSol_function(first,A,b);

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

function [slip] = multipleSol_function(stPt,A,b)


r = stPt.r;
B = stPt.B;
N = stPt.N;
D = stPt.D;
xb = stPt.xb;
yb = stPt.yb;
s = stPt.s;
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
        slip(jj).xb = xb1; slip(jj).yb = yb1; slip(jj).M = sum(xb1);slip(jj).s=s;
        
        jj = jj+1;
end
end
end
    
%%


    