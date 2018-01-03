%% Try1
% a = {[1 2], [1 2], [4 5]};
% allcomb(a{:})

% sets = {[1 2], [1 2], [4 5]};
% [x y z] = ndgrid(sets{:});
% cartProd = [x(:) y(:) z(:)];
clear;
close;

%% Try2 
% % sets = {[1 2], [1 2], [4 5]};

% % function cartesianProduct;
% 
% function result = cartesianProduct(sets)
%     c = cell(1, numel(sets));
%     [c{:}] = ndgrid( sets{:} );
%     result = cell2mat( cellfun(@(v)v(:), c, 'UniformOutput',false) );
% end
% % % cartProd = sortrows(cartProd, 1:numel(sets)); %you can sort the results:
% % % sets = cellfun(@unique, sets, 'UniformOutput',false);%check if the sets have no duplicate values
%% Try 3 combvec based method
% vectors = {[1 2], [3 6 9], [10 20]};
% combs = combvec(vectors{:});
% combs = sortrows(combvec(vectors{:}).')

%% Try 4 -NGrid based method is faster

% % % function combs = f1(vectors)
% % n = numel(vectors); %// number of vectors
% % combs = cell(1,n); %// pre-define to generate comma-separated list
% % [combs{end:-1:1}] = ndgrid(vectors{end:-1:1}); %// the reverse order in these two
% % %// comma-separated lists is needed to produce the rows of the result matrix in
% % %// lexicographical order
% % combs = cat(n+1, combs{:}); %// concat the n n-dim arrays along dimension n+1
% % combs = reshape(combs,[],n);

%% timing the above two codes--Check??

% nn = 20:20:240;
% t1 = [];
% t2 = [];
% for n = nn;
%     %//vectors = {1:n, 1:n, 1:n};
%     vectors = {1:n/10, 1:n, 1:n*10};
%     t = timeit(@() f1(vectors));
%     t1 = [t1; t];
%     t = timeit(@() f2(vectors));
%     t2 = [t2; t];
% end

%% try 5

% % % tic;
% % % v = {[1 2], [3 6 9], [10 20]};
% % % 
% % % L = [0 cumsum(cellfun(@length,v))];
% % % V = cell2mat(v);
% % % 
% % % J = nchoosek(1:L(end),length(v));
% % % J(any(J>repmat(L(2:end),[size(J,1) 1]),2) | ...
% % %   any(J<=repmat(L(1:end-1),[size(J,1) 1]),2),:)  = [];
% % % 
% % % V(J)
% % % toc

%% matfun nchoosek

v = 1:1:12;
C = nchoosek(v,5);
Cl = int2str(C);
B = [[1 2 3];[4 5 6];[7 8 9];[10 11 12]];
j =1;
k=1;

for i=1:1:792
    sum = 0;
    A = C(i,:);
     if ~(((numel(intersect(C(i,:),B(1,:)))  == 3))||(numel(intersect(C(i,:),B(2,:))) ==3)...
             ||(numel(intersect(C(i,:),B(3,:))) == 3)||(numel(intersect(C(i,:),B(4,:)))==3));
         
              C_new(j,:) = A;
%             
              if (numel(intersect(C_new(j,:),B(1,:))) == 1) sum=sum+1; end
              if (numel(intersect(C_new(j,:),B(2,:))) == 1)  sum=sum+1; end
              if (numel(intersect(C_new(j,:),B(3,:))) == 1)  sum=sum+1; end
              if (numel(intersect(C_new(j,:),B(4,:))) == 1)  sum=sum+1; end
              if sum ==3
                  C_f(k,:) = C_new(j,:);
                  k = k+1;
              end               
              j=j+1;
     end
end
    
 
 