function moments = CalcExactMoments(model_outputs)
% input is from running 
% 1. outputs = mexSolveModelGivenParms(parm_struc, useDevice, DisplayInterval)
% 2. model_outputs = ConvertCPPtoMAT(parm_struc, outputs)
    
    proby_cutoff = 1e-3;

    SPARSE_TRANSITION_PROB_YBD = GenSparseTransitionMatrix(p, model_outputs.defprob, model_outputs.idx_bnext, proby_cutoff);
    
    [V,D,FLAG] = eigs(SPARSE_TRANSITION_PROB_YBD',1);
    
    % checks
    if FLAG~=0 || abs(D(1,1)-1)>1e-8
        warning('moments:eigs may not have converged to the correct solution.');
    end
    
    prob_YBD = reshape(V(:,1)/sum(V(:,1)),[p.gridsize_tfp p.gridsize_b 2]);
    
    

end