function [gmatrix] = applySpin(spin,DCMtrx)

rotdcmat = spin2mat(spin); % matrix rotates -spin
gmatrix=  rotdcmat*DCMtrx;
%     dcg1=  rotdcmat*DCMtrx; %new axis rotation dc will be 
%     gmatrix = rotation('matrix',dcg1); 
end
    