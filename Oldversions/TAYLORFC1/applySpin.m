function [gmatrix] = applySpin(spin,DCMtrx)

rotdcmat = spin2mat(spin); % matrix rotates -spin ; rotdcmat is component tr matrix
    dcg1=  DCMtrx*rotdcmat'; %new axis rotation dc will be 
    gmatrix = rotation('matrix',dcg1); 
end
    