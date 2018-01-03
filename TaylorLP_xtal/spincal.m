
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