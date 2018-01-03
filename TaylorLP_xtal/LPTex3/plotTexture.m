load('OutputTexture.mat', 'new_gmatrixR')
cs = crystalSymmetry('-43m'); 
ss = specimenSymmetry('mmm');
o = orientation('Euler',new_gmatrixR(:,:,1), new_gmatrixR(:,:,2), new_gmatrixR(:,:,3),cs,ss);

%%
n = length(o);
                                        % n =20;
figure('position',[50 50 1200 500])
axesPos = subplot(5,5,1);
plotPDF(o(:,1),Miller(1,1,1,cs),'all',axesPos,'MarkerSize',3);
hold on
for i=1:1:n

axesPos = subplot(5,5,i);
plotPDF(o(:,i),Miller(1,1,1,cs),'all','parent',axesPos,'MarkerSize',3)
text(0,0,num2str(i-1));
end
hold off

