function [gmatrix] = applySpin(spin,DCMtrx)

rotdcmat = spin2mat(spin); % matrix rotates -spin
    dcg1=  DCMtrx*rotdcmat'; %new axis rotation dc will be 
    gmatrix = rotation('matrix',dcg1); 
end
    