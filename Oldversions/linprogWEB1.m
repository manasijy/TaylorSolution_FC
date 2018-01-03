%LINPROG  This code uses the revised simplex method to solve the linear
%       programming problem: Minimize the cost c'x subject to 
%       equations Ax=b and nonnegativity x >=0:
%
%                f =  Min { c'x; Ax = b, x >= 0 },
%
%       You must define an m by n matrix A , a column vector b with 
%       m components, and a column vector c with n components.
%%%       You may define k to change the first k equations of Ax=b to
%%%       inequalities: thus [A(1:k,:)]x <= b(1:k).   
%
%       The output vector x gives the minimum cost, which is the output f.
%
%                [x,f,itn] = linprog(A,b,c,k,maxit,tol)
%
%       At most "maxit" iterations (default 10*length(b)) are applied
%       and the actual number of iterations is returned in "itn".
%
%       If the optimal solution is unbounded or the constraints are
%       inconsistent then a diagnostic is displayed.  
%       Bland's rule is used to resolve degeneracies.  In exact
%       arithmetic cycling is not possible.  But in real life!!!
%       If x has more than 20 components it is returned in the
%       sparse format. Note that if A has many zeros it is worth
%       passing it to linprog in sparse format. 
%
%       Although written for teaching purposes this routine has
%       successfully solved some problems with size(A) = [50,100000]!
%
%       Please report any difficulties to: idc@math.canterbury.ac.nz

%       New version                  (c) I.D.Coope, 1988, 1993

%load(Aeq_Goss);
%[1;-0.5;1.0;0;0]

function [x,fval,it,B] = linprogWEB1(A,b,c,k,maxit,tol)
[m,n]=size(A); b=b(:); c=c(:); it=0;  %Converting b,c to column matrices
if (length(c)~=n || length(b)~=m),error('wrong dimensions'); end
if (nargin<6), tol=1e-10; end
if (nargin<5), maxit=10*m; end
if (nargin<4), k=0; elseif isempty(k), k=0; end
D=sign(sign(b)+.5); %This step is taken to drag out the signs of b.0.5 added to remove 0 
if k, D(1:k)=ones(k,1); end %For inequalities signs are taken as positive
D = diag(D);                    % initial (inverse) basis matrix
A = [A D];                      % incorporate slack/artificial variables.

%%
B = n+1:n+m;                    % initial basis it contains no.s of vars in basis eg x1,x4,x8 etc
N = 1:n;                        % non-basis it contains no.s of vars in nonbasis eg x2,x3,x5,x6,x7 etc
[bmin,j]=min(b(1:k));   % Why minimum of inequalities and its row is looked for???
if bmin<0,
   phase=1; xb=ones(m,1); s=[zeros(n+k,1);ones(m-k+1,1)];   % supercost
   % one element is added with 1. xb is of the size of no of eqns. s is a
   % col with (total no of vars + no of slacks) zeros followed by (no of
   % equalities +1)ones. i.e. one extra equality is added to the set
   N=[N,B(j)];  %Take the basis element corresponding to the -ve b inequality, and put that in non basis 
   J=B; J(j)=[]; %jth element of J(which is a copy of B) is removed.
   B(j)=n+m+1;  % and a new jth element of B (basis) 'n+m+1' is added 
   a=b-sum(A(:,J)')'; % a (Col Vect) = b (Col Vect) - sum(all (col vects)' of basis elements in J vect)'
   A=[A a]; % Now add one more column to A
   D(:,j)= -a/a(j);% jth coumn of D is -a/aj, this makes D(j,j)=-1
   D(j,j)=1/a(j); % jth element of Diag Mtx D is also updated. Why??
elseif k==m,        % i.e. all inequalities
   phase=2; xb=b; s=[c;zeros(m,1)];        % cost function
else                            % k< m and bmin>=0 -- This one is for Taylor
   phase=1; xb=abs(b); s=[zeros(n+k,1);ones(m-k,1)];   % supercost
end
while phase<3,
   df=-1; t=inf;
   yb= D'*s(B);  % multipliers for Ax=b. Took s(),(1 or 0), B'*pi  = c_B.It is using revised simplex pi = yb
   while (it < maxit)
      if isempty(N), break, end     % no freedom for minimization
      r = s(N) - [A(:,N)]'*yb;      % reduced costs: r_q = c_q - c_B'*inv(B)*N_q ;Or r_q = c_q-pi'*Nq
                                    %Here q are elements of N and ReducedCost is tested for all elements of N
      [rmin,q] = min(r);            % determine new basic variable % if reduced cost is zero,
                                    %there is no improvement i.e. another optimal solution
      if rmin>=-tol*(norm(s(N),inf)+1), break, end % optimal!
      it=it+1;
      if df>=0                      % apply Bland's rule to avoid cycling
         if maxit==inf,
            disp(['LINPROG(',int2str(it),'): warning! degenerate vertex']);
         end
         J=find(r<0); Nq=min(N(J)); q=find(N==Nq); %1-find the q for which r is -ve 2)-Nq is the first column with -ve r 
                                                   %3)- find the location of Nq in the array N          
      end 
      d = D*A(:,N(q)); % Descent Direction: B*d = -Nq.If dq >= 0, then the linear program is unbounded. Here d<=0 =>unbounded
      I=find(d>tol); %    I=find(d>0); Here it says that if d <0 then solution is bounded. This is because sign on RHS is changed above
      if isempty(I), disp('Solution is unbounded'); it=-it; break; end
      xbd=xb(I)./d(I); [r,p]=min(xbd); p=I(p);
      if df>=0,                     % apply Bland's rule to avoid cycling
         J=find(xbd==r); Bp=min(B(I(J))); p=find(B==Bp); 
      end 
      xb= xb - r*d; xb(p)=r;        % update x
      df=r*rmin;                    % change in f 
      v = D(p,:)/d(p);              % row vector
      yb= yb + v'*( s(N(q)) - d'*s(B) );
      d(p)=d(p)-1;
      D = D - d*v;                  % update inverse basis matrix
      t=B(p); B(p)=N(q);
      if t>n+k, N(q)=[]; else N(q)=t; end
   end                             % end of phase
   xb=xb+D*(b-A(:,B)*xb);          % iterative refinement
   I=find(xb<0);                   % must be due to rounding error
   if I, xb(I)=xb(I)-xb(I); end    % so correct
   if phase==2 | it<0, break; end; % B, xb,n,m,res=A(:,B)*xb-b
   if xb'*s(B)>tol,it=-it; disp('no feasible solution'); break;  end
   phase=phase+1;                  % re-initialise for Phase 2
   s=1e6*norm(c,'inf')*s; s(1:n)=c;% tol=tol*norm(s,inf);
end
x=sparse(n,1); x(B)=xb; x=x(1:n); if n<21, x=full(x); end
fval=c'*x;
if it>=maxit, disp('too many iterations'); it=-it; end




