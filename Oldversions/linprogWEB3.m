tic
clear

%%
gvector = [10,20,90];
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
% fprintf('Total Shear= %4.3f \n Shear: %c \n Iterations: %4.2d \n Basis: %d',fval,x,it,B);
% Output.x = x;
% Output.fval = fval;
% Output.it = it;
% Output.Basis = A(:,B);
% Output.B = B;
% Output.ReducedCost = r;
% Output.D = D;
toc


       %% Modification
       
zero_r = find(abs(r)<1e-3);
if ~isempty(zero_r)
% multi(1).x=sparse(n,1); multi(1).x(B)=xb; multi(1).x=x(1:n); if n<21, multi(1).x=full(multi(1).x); end       
slip(1).B = B; slip(1).shear = xb; slip(1).M = sum(xb);
jj = 1;
l_s = numel(zero_r);                                             
zero_d = D*A(:,N(zero_r));
round_d = round(zero_d);
all1 =  ones(5,1);
xzero = find(xb<0.001);
all1(xzero) = 0;


for i = 1:1:l_s
        
        mark = all1-round_d(:,i);
        if (min(mark) < 0), continue; end
        jj = jj+1;
        rp = find(round_d(:,i) == 1);
        [xmin,xminp] = min(xb(rp));
        % if (abs(xmin)>0.0001)
        B(rp(xminp)) = N(zero_r(i));
        mxb= xb - xmin*zero_d(:,i); mxb(rp(xminp))=xmin;
        slip(jj).B = B; slip(jj).shear = mxb; slip(jj).M = sum(mxb);
%         multi(i+1).x =sparse(n,1); multi(i+1).x(B) = mxb; multi(i+1).x=x(1:n); if n<21, multi(i+1).x=full(multi(i+1).x); end
end
end

    
 
    
    
    
      %%      
