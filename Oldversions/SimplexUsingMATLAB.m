% ebv is pivot column and lbv is the pivot row

function SimplexUsingMATLAB()
clc
format compact
format rational
A =[0.4 0.6 1 0 0 8.5;3 -1 0 1 0 25;3 6 0 0 1 70;-990 -900 0 0 0 5280];
[n m]= size(A);
fprintf('This is Table 1 - note the canonical form\n'),disp(A);
for iloop = 1:m
    [val ebv] = min(A(n,:));
    if val < 0
        fprintf('\nEntering Basic Variable Column  : '),disp(ebv);
        lbv = LBV(A,ebv);
        fprintf('\nLeaving Basic Variable (Pivot Row) : '),disp(lbv)
    else
        fprintf('\nNo further iterations possible')
        fprintf('\nYou may have the solution')
        fprintf('\nFinal table :')
        format short g;
        fprintf('\n-------------\n'),disp(A);
        break;
    end
    if lbv > 0
        A = RowOperations(A,ebv,lbv);
        textstr = strcat('Simplex Table :',num2str(iloop+1));
        disp(textstr);
        disp(A);
        
    end
end

function ret = LBV(A,i)
% i is the entering basic variable column
[n m] = size(A);
j = 0;  min = 1000;
for k = 1:n
    if A(k,i) > 0
        row = A(k,m)/A(k,i);
        if row < min
            min = row;  j = k;
        end
    end
end

if j == 0
    fprintf('Not possible to evaluate EBV ')
    RETURN
end
ret = j;

function ret = RowOperations(A,i,j)
% i is the EBV column
% j is the pivot row
[n m] = size(A);
A(j,:) = A(j,:)/A(j,i);
for k = 1: n
    if (k ~= j)
        A(k,:) = A(k,:) - A(k,i)*A(j,:);
    end
end
ret = A;
