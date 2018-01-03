function plot_plasticSpin(h,D) %(s,varargin)
% This function plots the poles of tensile axis on the sterogram as it is
% rotated due to shear on the given slip system 
%
% Input
%  D: an array of miller indices of rotated axis 
%
% Output


s = crystalSymmetry('-43m');
mtexFig = h;

%%
m = D;


% m = unique(m);

options = [{'MarkerSize',3,'MarkerFaceColor',[0.75 0 0.75],...
    'Marker','*','MarkerEdgeColor','k','grid','doNotDraw','antipodal'}];
  
hold(mtexFig.gca,'all')

for i = 1:length(m-1)
  m(i).scatter(options{:});
end

options = [{'MarkerSize',6,'MarkerFaceColor','g',...
    'Marker','d','MarkerEdgeColor','k','grid','doNotDraw','antipodal'}];
m(length(m)).scatter(options{:});




