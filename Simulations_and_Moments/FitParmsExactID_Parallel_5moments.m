function [model,stats,p,FittedParms,FitVal, EXITFLAG] = FitParmsExactID_Parallel_5moments(p_base, pool, ...
    target_debt_gdp, target_mean_spread, target_std_spread, target_recovery_frac, target_bid_ask_ND, ...
    beta_guess, frac_loss_y1_guess, slope_guess, recovery_bmax_guess, holding_cost_guess, ...
    DisplayInterval,T,TBurn,Us_z,Us_m,Us_reenter,QuarterLength,NumQuartersExcludeAfterReentry)

% use same shocks to get rid of simulation noise
% moments are fixed: debt level, mean spreads, vol of spreads.
% parameters: beta, level of costs, slope of costs

    NumGPDDevices = gpuDeviceCount;

    if pool.NumWorkers>NumGPDDevices
        error('Insufficient number of GPU devices');
    end
    
    % DiffMinChange: Minimum change in variables for finite-difference gradients (a positive scalar). The default is 0.
    % TypicalX
    
    TypX = ones(5,1);
    % fsolve_opts = optimset('display','iter-detailed','TypicalX',TypX,'DiffMinChange',0.01,'TolFun',1e-2,'UseParallel',true);
    fsolve_opts = optimset('display','iter-detailed','TypicalX',TypX,'DiffMinChange',0.01,'TolFun',1e-2,'TolX',1e-5,'UseParallel',true);
    
    if beta_guess<=0.97; error('check beta_guess'); end
    
    X0 = zeros(5,1);
    X0(1) = -log(0.03/(beta_guess - 0.97) - 1)/5 + 1;
    X0(2) = frac_loss_y1_guess/0.07;
    X0(3) = slope_guess/0.3;
    X0(4) = 1 + log(recovery_bmax_guess/0.8);
    X0(5) = 1 + log(holding_cost_guess/0.005);

    [XFit, FitVal, EXITFLAG] = fsolve(@(X)Moments(X),X0,fsolve_opts);
    
    loss_y1_fit = 0.07*XFit(2); % typical value around 0.07
    loss_slope_fit = 0.3*XFit(3); % typical value around 0.3
    
    FittedParms.beta = 0.97 + 0.03/(1 + exp(5*(1 - XFit(1))));
    FittedParms.d_y = loss_y1_fit - loss_slope_fit;
    FittedParms.d_yy = loss_slope_fit;
    FittedParms.recovery_bmax = 0.8*exp(XFit(4) - 1);
    FittedParms.holding_cost = 0.005*exp(XFit(5)-1);

    function Z = Moments(X)
        
        t = getCurrentTask();
        if isempty(t)
            useDevice = 0;
        else
            if t.ID<=NumGPDDevices
                useDevice = t.ID - 1;
            else
                fprintf('t.ID=%g exceeds NumGPDDevices of %g.\n', t.ID, NumGPDDevices);
                error('t.ID exceeds number of workers');
            end
        end
        
        Z = zeros(5,1);
        p = p_base;
        
        % assign parameters:
        % scale input parameters so that they have magnitude 1 and a
        % reasonable value for DiffMinChange is 0.01

        loss_y1 = 0.07*X(2); % typical value around 0.07
        loss_slope = 0.3*X(3); % typical value around 0.3
        
        p.beta = 0.97 + 0.03/(1 + exp(5*(1 - X(1)))); % typical value around 0.985
        p.d_y = loss_y1 - loss_slope;
        p.d_yy = loss_slope;
        p.recovery_bmax = 0.8*exp(X(4) - 1);
        p.holding_cost = 0.005*exp(X(5)-1);
        
        fprintf('beta=%g, d_y=%g, d_yy=%g, recovery_bmax=%g, holding_cost=%g.\n', p.beta, p.d_y, p.d_yy, p.recovery_bmax, p.holding_cost);
        
        % solve model on GPU
        model_mex = mexSolveModelGivenParms(p, useDevice, DisplayInterval);
        
        % convert from c++ output to matlab output
        model = ConvertCPPtoMAT(p, model_mex);
        
        % simulate
        [~,stats] = SimSinglePath_mex_output(p,model,T,TBurn,Us_z,Us_m,Us_reenter,QuarterLength,NumQuartersExcludeAfterReentry);
        
        % output, assumes TolFun is 1e-2
        Z(1) = 10*(stats.DebtToOutput_mean - target_debt_gdp)/target_debt_gdp;
        Z(2) = 10*(stats.Spread_mean - target_mean_spread)/target_mean_spread;
        Z(3) = 10*(stats.Spread_std - target_std_spread)/target_std_spread;
        Z(4) = 10*(stats.RecoveryFraction_mean - target_recovery_frac)/target_recovery_frac;
        Z(5) = 10*(stats.bidask_ND_mean - target_bid_ask_ND)/target_bid_ask_ND;
        
        fprintf('Moment(1)=%g, Moment(2)=%g, Moment(3)=%g, Moment(4)=%g, Moment(5)=%g.\n', ...
            stats.DebtToOutput_mean, stats.Spread_mean, stats.Spread_std, stats.RecoveryFraction_mean, stats.bidask_ND_mean);
        
    end

end