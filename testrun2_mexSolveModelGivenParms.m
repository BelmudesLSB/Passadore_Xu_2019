clear

%% read

load('../../../liquidity.nov.2017.tauchen.maxdefprob/calibration_final/BASELINE.mat','p','model')

p.prob_meet_mkt_maker_ND = p.prob_meet_mkt_maker;
p.prob_meet_mkt_maker_D = p.prob_meet_mkt_maker;
p.tfp_lb = log(p.grid_y_nd(1));
p.tfp_ub = log(p.grid_y_nd(end));
p.conv_chk_tfp_lb = p.tfp_lb;
p.conv_chk_tfp_ub = p.tfp_ub;
p.conv_chk_b_ub = p.bmax;

useDevice = 0;
DisplayInterval = 100;

model_mex = mexSolveModelGivenParms(p, useDevice, DisplayInterval);

%% compare solutions
% mex solution is stored in c++ format over (idxy*GRIDSIZE_B + idxb)
% previous solution is in matlab's format over (y,b)

err_CVD = max(abs(reshape(permute(reshape(model_mex.CVD,[p.gridsize_b p.gridsize_tfp]),[2 1]),[],1) - model.CVD(:)))
err_CVND = max(abs(reshape(permute(reshape(model_mex.CVND,[p.gridsize_b p.gridsize_tfp]),[2 1]),[],1) - model.CVND(:)))
err_QHND = max(abs(reshape(permute(reshape(model_mex.qH_ND,[p.gridsize_b p.gridsize_tfp]),[2 1]),[],1) - model.qNDH(:)))
err_QLND = max(abs(reshape(permute(reshape(model_mex.qL_ND,[p.gridsize_b p.gridsize_tfp]),[2 1]),[],1) - model.qNDL(:)))
err_QHD = max(abs(reshape(permute(reshape(model_mex.qH_D,[p.gridsize_b p.gridsize_tfp]),[2 1]),[],1) - model.qDH(:)))
err_QHL = max(abs(reshape(permute(reshape(model_mex.qL_D,[p.gridsize_b p.gridsize_tfp]),[2 1]),[],1) - model.qDL(:)))
err_defprob = max(abs(reshape(permute(reshape(model_mex.defprob,[p.gridsize_b p.gridsize_tfp]),[2 1]),[],1) - model.defprob(:)))
err_defthresh = max(abs(reshape(permute(reshape(model_mex.defthresh,[p.gridsize_b p.gridsize_tfp]),[2 1]),[],1) - model.defthresh(:)))

idx_region = cellfun(@(X)find(X.eps_mesh(1:end-1)<=p.ggq_mmean & p.ggq_mmean<=X.eps_mesh(2:end),1,'first'),model.policy);
idx_bnext_old = arrayfun(@(i)model.policy{i}.idx_bnext(idx_region(i)),1:numel(model.policy));
err_idxbnext = max(abs(reshape(permute(reshape(model_mex.idx_bchoice+1,[p.gridsize_b p.gridsize_tfp]),[2 1]),[],1) - idx_bnext_old(:)))
