function [gmatrix] = applySpin(spin,DCMtrx)

rotdcmat = spin2mat(spin);
    dcg1=  rotdcmat'*DCMtrx'; 
    gmatrix = rotation('matrix',dcg1); 
end
    