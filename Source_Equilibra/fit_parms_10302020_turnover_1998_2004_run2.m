clear all

%% base parameters

p.rra = 2;
p.beta = 0.9842;
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
p.holding_cost = 0.0058;
p.bargain_power_mkt_maker_ND = 0.5;
p.bargain_power_mkt_maker_D = 0.5;
p.reenterprob = 0.0128;
p.debt_z = 0.01;
p.debt_m = 0.0167;
p.debt_writedown_reentry = 0;
p.max_def_prob = 0.75;
p.recovery_bmax = 0.83;
p.gridsize_tfp = 200;
p.gridsize_b = 450;
p.bmax = 6;
p.grid_b = linspace(0,p.bmax,p.gridsize_b)';
p.rho_tfp = 0.9830;
p.sigma_tfp = 0.0151;

TURNOVER_MONTHLY = 1.19/12; % average: emta volume sovereign eurobonds (non-brady) divided by ppg.

p.prob_meet_mkt_maker_ND = 0.8647; % = 1 - exp(-4). So average trade is 1 week
p.prob_meet_mkt_maker_D = 0.8647;

prob_liqshock = TURNOVER_MONTHLY*(p.debt_m+(1-p.debt_m)*p.prob_meet_mkt_maker_ND)/(p.prob_meet_mkt_maker_ND-TURNOVER_MONTHLY*(1-p.debt_m));

if prob_liqshock<=0; error('prob_liqshock<=0.'); end

p.prob_liqshock = prob_liqshock;
uncond_std = p.sigma_tfp/sqrt(1 - p.rho_tfp^2);
width_numstd = 3;
width_numstd_conv_chk = 2.5;
p.tfp_lb = -width_numstd*uncond_std;
p.tfp_ub = width_numstd*uncond_std;
p.conv_chk_tfp_lb = -width_numstd_conv_chk*uncond_std;
p.conv_chk_tfp_ub = width_numstd_conv_chk*uncond_std;
p.conv_chk_b_ub = 5.5;
p.d_y = -0.2530;
p.d_yy = 0.3241;
p.d_b = 0;
p.d_bb = 0;
p.d_yb = 0;
p.maxiters = 8000;
p.errtolV = 1e-6;
p.errtolQ = 1e-5;
p.NumUpdateRegimes = 7;
p.UpdateRegimeIter = [200;500;1000;2000;4000;6000];
p.UpdateWeightOldV = zeros(7,1);
p.UpdateWeightOldQ = [0;0.5;0.9;0.95;0.98;0.98;0.98];
p.UpdateWeightOldDefProb = [0;0;0;0;0.9;0.95;0.98];

%% specify moments

target_debt_gdp = 1.0;
target_mean_spread = 0.0815;
target_std_spread = 0.0443;
target_recovery_frac = 0.3;
target_bid_ask_ND = 0.0056;

beta_guess = 0.984239998133106;
d_y = -0.261760012535829;
d_yy = 0.334781605794710; 
frac_loss_y1_guess = d_y + d_yy;
slope_guess = d_yy;
recovery_bmax_guess = 0.809952250437929; 
holding_cost_guess = 0.005503259285985;

%% run fitting procedure

% Specs of the GPU:
poolsize = gpuDeviceCount;
pool = parpool(poolsize, 'IdleTimeout', 240);

DisplayInterval = -1; % don't display
TBurn = 12*10^3;
T = TBurn + 12*10^6;

rng(0); % set seed to ensure reproducibility
Us_z = rand(T,1);
Us_m = rand(T,1);
Us_reenter = rand(T,1);
QuarterLength = 3;
NumQuartersExcludeAfterReentry = 20;

tic
[model,stats,p_fit,FittedParms,FitVal, EXITFLAG] = FitParmsExactID_Parallel_5moments(p, pool, ...
    target_debt_gdp, target_mean_spread, target_std_spread, target_recovery_frac, target_bid_ask_ND, ...
    beta_guess, frac_loss_y1_guess, slope_guess, recovery_bmax_guess, holding_cost_guess, ...
    DisplayInterval,T,TBurn,Us_z,Us_m,Us_reenter,QuarterLength,NumQuartersExcludeAfterReentry);
toc

delete(pool)
clear Us_z Us_m Us_reenter

%% save

save([mfilename,'.mat'])