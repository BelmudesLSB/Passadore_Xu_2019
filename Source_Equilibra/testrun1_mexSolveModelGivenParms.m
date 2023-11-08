clear

%% specify parameters

p.rra = 2;
p.beta = 0.9841;
p.clb = 1e-6;

p.ggq_mstd = 0.004;
p.ggq_mlb = -3*p.ggq_mstd;
p.ggq_mub = 3*p.ggq_mstd;
p.ggq_mmean = 0.5*(p.ggq_mub + p.ggq_mlb);
p.ggq_mDval = p.ggq_mmean;

p.errtol_bisect = 1e-12;
p.maxiters_bisect = 1000;
p.rH = 0.0033;
p.rL = 0.0033;
p.holding_cost = 0.0014;
p.prob_liqshock = 0.1393;
p.prob_meet_mkt_maker_ND = 0.8647;
p.prob_meet_mkt_maker_D = 0.8647;
p.bargain_power_mkt_maker_ND = 0.8750;
p.bargain_power_mkt_maker_D = 1;
p.reenterprob = 0.0128;
p.debt_z = 0.01;
p.debt_m = 0.0167;
p.debt_writedown_reentry = 0;
p.max_def_prob = 0.75;
p.recovery_bmax = 0.83;
p.gridsize_tfp = 32;
p.gridsize_b = 128;
p.bmax = 6;
p.rho_tfp = 0.9830;
p.sigma_tfp = 0.0151;

uncond_std = p.sigma_tfp/sqrt(1 - p.rho_tfp^2);
width_numstd = 3;
width_numstd_conv_chk = 2.5;

p.tfp_lb = -width_numstd*uncond_std;
p.tfp_ub = width_numstd*uncond_std;
p.conv_chk_tfp_lb = -width_numstd_conv_chk*uncond_std;
p.conv_chk_tfp_ub = width_numstd_conv_chk*uncond_std;
p.conv_chk_b_ub = 5.5;

p.d_y = -0.2640;
p.d_yy = 0.3370;
p.d_b = 0;
p.d_bb = 0;
p.d_yb = 0;
p.maxiters = 3000;
p.errtolV = 1e-5;
p.errtolQ = 1e-5;
p.NumUpdateRegimes = 7;
p.UpdateRegimeIter = [200;500;1000;2000;4000;6000];
p.UpdateWeightOldV = zeros(7,1);
p.UpdateWeightOldQ = [0;0.5;0.9;0.95;0.98;0.98;0.98];
p.UpdateWeightOldDefProb = [0;0;0;0;0.9;0.95;0.98];

%% read

useDevice = 0;
DisplayInterval = 100; % <1: don't display
% DisplayInterval = -1; % <1: don't display

tic
p_out = mexSolveModelGivenParms(p, useDevice, DisplayInterval)
toc