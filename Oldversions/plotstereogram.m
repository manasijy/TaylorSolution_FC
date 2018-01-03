function [h] = plotstereogram() %(s,varargin)
% This program plots the standard cubic (001) sterogram



s = crystalSymmetry('-43m');
mtexFig = newMtexFigure('xAxisDirection','east'); %(varargin{:});

%% which directions to plot
m = [Miller(1,0,0,s),Miller(0,1,0,s),...
  Miller(0,0,1,s),Miller(1,1,0,s),...
  Miller(0,1,1,s),Miller(1,0,1,s),...
  Miller(1,1,1,s)];


m = unique(m);
options = [{'symmetrised','labeled','MarkerEdgeColor','k','grid','doNotDraw','antipodal','xAxisDirection','south'}]; %
  
% plot them

washold = getHoldState(mtexFig.gca);
hold(mtexFig.gca,'all')
for i = 1:length(m)
  m(i).scatter(options{:});
end
hold(mtexFig.gca,washold)


% postprocess figure
setappdata(gcf,'CS',s);
set(gcf,'tag','ipdf');
mtexFig.drawNow('figSize',getMTEXpref('figSize'));
h = mtexFig;
