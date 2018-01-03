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
    
    
    