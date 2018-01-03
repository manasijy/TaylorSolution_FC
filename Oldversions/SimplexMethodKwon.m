function [xsol objval exitflag]=SimplexMethodKwon(c, Aeq, beq, B_set)
% Simplex_Method solves a linear program in standard form:
% min c'*x
% s.t. Aeq*x = beq
% x >= 0
% by using the simplex method of George B. Dantzig
%
% Inputs:
% c = n*1 vector of objective coefficients
% Aeq = m*n matrix with m < n
% beq = m*1 vector of right hand side (RHS) coefficients
% B_set = m*1 vector that contains indices (subscripts) of basic variables
%
% Parameter and Variable Partitions
%
% c, Aeq, and beq are partitioned according to partition of x into
% x' =[x_B' | x_N'] where
% x_B is an m*1 vector of basic variables
% x_N is an (n-m)*1 vector of non-basic variables
% c' = [c_B' | c_N']
% c_B is the objective coefficients of x_B, an m*1 vector
% c_N is the objective coefficients of x_N, an (n-m)*1 vector
% Aeq = [B | N]
% B is the m*m basis matrix
% N is m*(n-m) non-basis matrix
% set = [B_set' | N_set']
% set is a set of indices (subscripts) of x
% N_set is an (n-m)*1 vector of indices (subscripts) of non-basic variables
%
% Output:
% xsol = n*1 vector, contains final solution of LP
% objval is a scalar, final objective value of LP
% iter is a struct that includes for every iteration the following details:
% B_Set, N_set, c_B, x_B, r_N, step, and d where step is a step length
% and d is a search direction.
%
% exitflag describes the exit condition of the problem as follows:
% 0 - optimal solution
% 1 - unbounded problem
xsol=[]; objval=[]; exitflag=[];
%% Step 0: Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate an initial basic feasible solution and partition c, x, and Aeq
% so that c=[c_B | c_N] x=[x_B | x_N] Aeq=[B | N]

set=[1:length(c)]';
set(find(ismember(set, B_set)==1))=[];
N_set=set;
B=Aeq(:,B_set); %basis matrix B
c_B=c(B_set); %obj coefficients of current basic variables
x_B=B\beq; %compute basic variables
N=Aeq(:,N_set); %non-basis matrix N
c_N=c(N_set); %obj coefficients of current non-basic variables
x_N=zeros(length(N_set),1); %x_N, non-basic variables equal 0
x=[x_B; x_N]; %partition x according to basis
obj=[c_B; c_N]'*x; %initial objective function value
k=0;
while k>=0
%% Step 1: Optimality Check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the reduced costs r_q=c_q-c_B'*B^(-1)*N_q for q in N_set
% if r_q >= 0, STOP, current solution optimal, else go to STEP 2
pie=B'\c_B; %solve the system B^T*pie=c_B for simplex multipliers
r_N=c_N'-pie'*N; % compute reduced cost for non-basic variables
ratioflag=find(r_N<0);
if isempty(ratioflag) %if r_q >= 0, then STOP. Optimal
disp('probelm solved')
exitflag=0;
objval=obj;
%subscripts of x are in ascending order
set_temp=[B_set; N_set];
for a=1:length(c)
xsol(a,1)=x(find(set_temp==a));
end
break
else % if r_q < 0, GO TO Step 2
%% Step 2: Descent Direction Generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct d_q=[-B^(-1)*N; e_q].
% If d_q >= 0, then LP is unbounded, STOP, else go to STEP 3.
enter=ratioflag(1); %choosing entering variable
e=zeros(length(N_set),1);
e(enter)=1; %construct vector e_q
d=-B\N(:,enter); %solve the system Bd=-N_q
direction=[d; e]; %improved direction d
d_flag=find(direction < 0);
if isempty(d_flag) %if direction > 0, then STOP.(unbounded)
disp('unbounded problem')
exitflag=1;
break

else %if d_q < 0, GO TO Step 3
%% Step 3: Step Length Generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute step length by the minimum ratio test. Go to STEP 4.
step_set=-x(d_flag)./direction(d_flag);
step=min(step_set);
%% Step 4: Improved Solution Generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Let x(k+1) = x(k) + alpha*d_q. Go to Step 5.
x_d=x+step*direction;
leave_set=find(x_d(1:length(B_set))==0);
leave=leave_set(1); %determining leaving variable
%% Step 5: Basis Update
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the new basis B for next iteration,
% Update c=[c_B|c_N],x=[x_B|x_N],& Aeq=[B|N]. Go to STEP 1.
B_set_temp=B_set;
N_set_temp=N_set;
x_B=x_d(1:length(B_set));
x_B_temp=x_d(1:length(B_set));
x_N_temp=x_d(length(B_set)+1:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%exchange the entering and leaving variables in B_set
B_set(find(B_set_temp==B_set_temp(leave)))=N_set_temp(enter);
N_set(find(N_set_temp==N_set_temp(enter)))=B_set_temp(leave);
x_B(find(x_B_temp==x_B_temp(leave)))=x_N_temp(enter);
B=Aeq(:,B_set); %update basis B
c_B=c(B_set); %update c_B
N=Aeq(:,N_set); %update non-basis N
c_N=c(N_set); %update c_N
x=[x_B; x_N]; %update x = [x_B | x_N]
obj=[c_B; c_N]'*x; %new objective value
k=k+1; %GO TO Step 1
end
end
end