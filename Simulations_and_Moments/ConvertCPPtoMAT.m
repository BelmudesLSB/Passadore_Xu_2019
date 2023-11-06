function conv_outputs = ConvertCPPtoMAT(p,model_outputs)
% input is from running 
% outputs = mexSolveModelGivenParms(parm_struc, useDevice, DisplayInterval)
    
    conv_outputs = model_outputs;
    conv_outputs.CVD = permute(reshape(model_outputs.CVD,[p.gridsize_b p.gridsize_tfp]),[2 1]);
    conv_outputs.CVND = permute(reshape(model_outputs.CVND,[p.gridsize_b p.gridsize_tfp]),[2 1]);
    conv_outputs.qNDH = permute(reshape(model_outputs.qH_ND,[p.gridsize_b p.gridsize_tfp]),[2 1]);
    conv_outputs.qNDL = permute(reshape(model_outputs.qL_ND,[p.gridsize_b p.gridsize_tfp]),[2 1]);
    conv_outputs.qDH = permute(reshape(model_outputs.qH_D,[p.gridsize_b p.gridsize_tfp]),[2 1]);
    conv_outputs.qDL = permute(reshape(model_outputs.qL_D,[p.gridsize_b p.gridsize_tfp]),[2 1]);
    conv_outputs.defprob = permute(reshape(model_outputs.defprob,[p.gridsize_b p.gridsize_tfp]),[2 1]);
    conv_outputs.defthresh = permute(reshape(model_outputs.defthresh,[p.gridsize_b p.gridsize_tfp]),[2 1]);
    conv_outputs.idx_bchoice = permute(reshape(model_outputs.idx_bchoice+1,[p.gridsize_b p.gridsize_tfp]),[2 1]);
    conv_outputs.prob_y = permute(reshape(model_outputs.prob_y,[p.gridsize_tfp p.gridsize_tfp]),[2 1]);
    conv_outputs.grid_y_d = permute(reshape(model_outputs.grid_y_d,[p.gridsize_b p.gridsize_tfp]),[2 1]);

end