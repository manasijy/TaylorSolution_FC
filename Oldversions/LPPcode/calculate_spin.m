function [spin,UniqSolution] = calculate_spin(Solution,criteria)

load('SlipSystem24.mat','SlipSystem');

        switch criteria 
            case 'minvar'
                UniqSolution = uniqueSol_minVar(Solution);
            case 'maxvar'
                UniqSolution = uniqueSol_maxVar(Solution);
            case 'minplasticspin'
                UniqSolution = uniqueSol_MinPlasticSpin(Solution);
        end
    SS= SlipSystem(UniqSolution.B);
    shear = UniqSolution.xb;
    spin = zeros(3,3);

    for ii=1:1:5, 
        spin = spin + shear(ii)*SS(ii).q.M;
    end