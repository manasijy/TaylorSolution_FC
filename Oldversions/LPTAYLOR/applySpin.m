function [gmatrix] = applySpin(spin_min,DCMtrx)

rotdcmat = spin2mat(spin_min);
    dcg1=  rotdcmat'*DCMtrx'; 
    gmatrix = rotation('matrix',dcg1); 
end
    