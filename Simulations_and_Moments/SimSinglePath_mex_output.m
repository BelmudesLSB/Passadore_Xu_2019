function [sim,stats] = SimSinglePath_mex_output(p,model,T,TBurn,Us_z,Us_m,Us_reenter,QuarterLength,NumQuartersExcludeAfterReentry)

    % one period ahead default prob
    if nargin<8
        QuarterLength = 3;
    end
    
    if nargin<9
        NumQuartersExcludeAfterReentry = 20;
    end

    if ~isfield(p,'cumul_prob_y')
        p.cumul_prob_y = cumsum(model.prob_y')';
    end
    
    if ~isfield(p,'cumul_stat_prob_y')
        [V,D] = eig(model.prob_y');
        [~,idx] = min(abs(diag(D)-1));
        stat_prob_y = V(:,idx)/sum(V(:,idx));
        p.cumul_stat_prob_y = cumsum(stat_prob_y);
    end
    
    T = QuarterLength*floor(T/QuarterLength);
    TBurn = QuarterLength*floor(TBurn/QuarterLength);

    yt = NaN*zeros(T,1);
    bt = NaN*zeros(T,1);
    Accesst = NaN*zeros(T,1); % Credit access indicator, prior to default decision, but after reentry shocks are drawn
    Dt = NaN*zeros(T,1); % D(t) = 0,1,NaN
    qHt = NaN*zeros(T,1);
    qLt = NaN*zeros(T,1);
    ct = NaN*zeros(T,1); % consumption
    NXt = NaN*zeros(T,1); % net exports
    bidaskt = NaN*zeros(T,1); % bid ask
    fracHt = NaN*zeros(T,1); % fraction bond holders of H type, start of period
    issuet = NaN*zeros(T,1);
    issuefract = NaN*zeros(T,1);
    defprob1t = NaN*zeros(T,1);
    Autarkyt = NaN*zeros(T,1);
    recovery_frac = NaN*zeros(T,1);
    
    X_idx_y = find(Us_z(1)<=p.cumul_stat_prob_y,1,'first');
    X_idx_b = 1;
    X_C1 = true;
    X_fracH = 1; % initially, all bonds issued to H types
    
    Plb = normcdf(p.ggq_mlb,p.ggq_mmean,p.ggq_mstd);
    Pub = normcdf(p.ggq_mub,p.ggq_mmean,p.ggq_mstd);
    et = norminv(Plb+Us_m*(Pub-Plb),p.ggq_mmean,p.ggq_mstd);
    
    for ixt = 1:T
        fracHt(ixt) = X_fracH;
        % reentry
        if ~X_C1
            if Us_reenter(ixt)<=p.reenterprob
                C2 = true;
                bdef = p.grid_b(X_idx_b);
                [~,X_idx_b] = min(abs(p.grid_b - min((1-p.debt_writedown_reentry)*p.grid_b(X_idx_b),p.recovery_bmax)));
                if bdef>0
                    recovery_frac(ixt) = p.grid_b(X_idx_b)/bdef;
                end
            else
                C2 = false;
            end
        else
            C2 = true;
        end
        
        bt(ixt) = p.grid_b(X_idx_b);
        X_idx_y = find(Us_z(ixt)<=p.cumul_prob_y(X_idx_y,:),1,'first');
        Accesst(ixt) = C2;
        
        if C2
            % check for default
            % idx_bnext=0: empty choice set
            if et(ixt)<=model.defthresh(X_idx_y,X_idx_b) || model.idx_bchoice(X_idx_y,X_idx_b)<=0
                bAutarky = true;
                Dt(ixt) = true;
            else
                bAutarky = false;
                Dt(ixt) = false;
            end
        else
            et(ixt) = p.ggq_mDval;
            bAutarky = true;
        end
        
        Autarkyt(ixt) = bAutarky;
        % policies
        if bAutarky
            X_C1 = false;
            yt(ixt) = model.grid_y_d(X_idx_y,X_idx_b) + et(ixt);
            qHt(ixt) = model.qDH(X_idx_y,X_idx_b);
            qLt(ixt) = model.qDL(X_idx_y,X_idx_b);
            ct(ixt) = yt(ixt);
            NXt(ixt) = 0;
            qsale = p.bargain_power_mkt_maker_D*qLt(ixt) + (1-p.bargain_power_mkt_maker_D)*qHt(ixt);
            bidaskt(ixt) = 2*(qHt(ixt) - qsale)/(qsale + qHt(ixt));
            X_fracH = (1-p.prob_liqshock)*X_fracH + p.prob_meet_mkt_maker_D*(1 - X_fracH);
        else
            X_C1 = true;
            yt(ixt) = model.grid_y_nd(X_idx_y) + et(ixt);
            X_idx_bnext = model.idx_bchoice(X_idx_y,X_idx_b);
            qHt(ixt) = model.qNDH(X_idx_y,X_idx_bnext);
            qLt(ixt) = model.qNDL(X_idx_y,X_idx_bnext);
            defprob1t(ixt) = model.defprob(X_idx_y,X_idx_bnext);
            issuet(ixt) = p.grid_b(X_idx_bnext) - (1-p.debt_m)*p.grid_b(X_idx_b);
            issuefract(ixt) = issuet(ixt)/(p.grid_b(X_idx_b));
            NXt(ixt) = (p.debt_m + (1-p.debt_m)*p.debt_z)*p.grid_b(X_idx_b) - issuet(ixt)*qHt(ixt);
            if p.grid_b(X_idx_bnext)>(1-p.debt_m)*p.grid_b(X_idx_b)
                X_fracH = (((1-p.prob_liqshock)*X_fracH + p.prob_meet_mkt_maker_ND*(1 - X_fracH))*(1-p.debt_m)*p.grid_b(X_idx_b) ...
                    + (1-p.prob_liqshock)*(p.grid_b(X_idx_bnext)-(1-p.debt_m)*p.grid_b(X_idx_b)))/p.grid_b(X_idx_bnext);
            else
                if X_idx_bnext==1
                    % bought everything back
                    X_fracH = 1;
                else
                    X_fracH = (1-p.prob_liqshock)*X_fracH + p.prob_meet_mkt_maker_ND*(1 - X_fracH);
                end
            end
            ct(ixt) = yt(ixt) - NXt(ixt);
            X_idx_b = X_idx_bnext;
            qsale = p.bargain_power_mkt_maker_ND*qLt(ixt) + (1-p.bargain_power_mkt_maker_ND)*qHt(ixt);
            bidaskt(ixt) = 2*(qHt(ixt) - qsale)/(qsale + qHt(ixt)); 
        end 
    end
    
    yt = yt(TBurn+1:end);
    bt = bt(TBurn+1:end);
    Accesst = Accesst(TBurn+1:end);
    Autarkyt = Autarkyt(TBurn+1:end);
    recovery_frac = recovery_frac(TBurn+1:end);
    Dt = Dt(TBurn+1:end); % D(t) = 0,1,NaN
    qHt = qHt(TBurn+1:end);
    qLt = qLt(TBurn+1:end);
    defprob1t = defprob1t(TBurn+1:end);
    ct = ct(TBurn+1:end); % consumption
    NXt = NXt(TBurn+1:end); % net exports
    bidaskt = bidaskt(TBurn+1:end); % bid ask
    fracHt = fracHt(TBurn+1:end);
    issuet = issuet(TBurn+1:end);
    issuefract = issuefract(TBurn+1:end);
    spreadHt = (1 + (p.debt_m + (1-p.debt_m)*p.debt_z)./qHt - p.debt_m).^(4*QuarterLength) ...
        - (1+p.rH)^(4*QuarterLength);
    spreadLt = (1 + (p.debt_m + (1-p.debt_m)*p.debt_z)./qLt - p.debt_m).^(4*QuarterLength) ...
        - (1+p.rH)^(4*QuarterLength);
    yt = reshape(sum(reshape(yt,QuarterLength,[]),1),[],1);
    bt_sum = reshape(sum(reshape(bt,QuarterLength,[]),1),[],1);
    bt = bt(QuarterLength:QuarterLength:end); % quarter end
    Accesst = Accesst(QuarterLength:QuarterLength:end); % quarter end
    Autarkyt = Autarkyt(QuarterLength:QuarterLength:end); % quarter end
    Dt = reshape(any(reshape(Dt,QuarterLength,[]),1),[],1); % any time during quarter
    qHt = reshape(mean(reshape(qHt,QuarterLength,[]),1),[],1); % average within quarter
    qLt = reshape(mean(reshape(qLt,QuarterLength,[]),1),[],1); % average within quarter
    spreadHt = reshape(mean(reshape(spreadHt,QuarterLength,[]),1),[],1); % average within quarter
    spreadLt = reshape(mean(reshape(spreadLt,QuarterLength,[]),1),[],1); % average within quarter
    defprob1t = reshape(mean(reshape(defprob1t,QuarterLength,[]),1),[],1); % average within quarter
    bidaskt = reshape(mean(reshape(bidaskt,QuarterLength,[]),1),[],1); % average within quarter
    fracHt = reshape(mean(reshape(fracHt,QuarterLength,[]),1),[],1); % average within quarter
    ct = reshape(sum(reshape(ct,QuarterLength,[]),1),[],1);
    NXt = reshape(sum(reshape(NXt,QuarterLength,[]),1),[],1);
    issuet = reshape(sum(reshape(issuet,QuarterLength,[]),1),[],1);
    issuefract = reshape(mean(reshape(issuefract,QuarterLength,[]),1),[],1);
    
    % store
    sim.yt = yt;
    sim.bt = bt;
    sim.Accesst = Accesst;
    sim.Autarkyt = Autarkyt;
    sim.Dt = Dt;
    sim.qHt = qHt;
    sim.qLt = qLt;
    sim.spreadHt = spreadHt;
    sim.spreadLt = spreadLt;
    sim.defprob1t = defprob1t;
    sim.ct = ct;
    sim.NXt = NXt;
    sim.bidaskt = bidaskt;
    sim.fracHt = fracHt;
    sim.issuet = issuet;
    sim.issuefract = issuefract;

    if NumQuartersExcludeAfterReentry==0
        idx_ND = find(~Autarkyt); % has access and didn't default
        stats.DefFreq = 4*nansum(Dt)/nansum(Accesst);% defined as the number of defaults per year.
    else
        bAdmit = zeros(size(Autarkyt));
        bAdmit(:) = false;
        for i = 1+NumQuartersExcludeAfterReentry:numel(Autarkyt)
            if all(~Autarkyt(i-NumQuartersExcludeAfterReentry:i))
                bAdmit(i) = true;
            end
        end
        idx_ND = find(bAdmit);
        
        NumDefs = 0;
        NumAllowable = 0;
        for i = 1+NumQuartersExcludeAfterReentry:numel(Autarkyt)
            if all(Accesst(i-NumQuartersExcludeAfterReentry:i))
                NumAllowable = NumAllowable + 1;
                if Dt(i)==true
                    NumDefs = NumDefs + 1;
                end
            end
        end
        stats.DefFreq = 4*NumDefs/NumAllowable;% defined as the number of defaults per year.
    end
    
    idx_D = find(Autarkyt);
    stats.DebtToOutput_mean = mean(bt(idx_ND)./yt(idx_ND));
    spread = (1 + (p.debt_m + (1-p.debt_m)*p.debt_z)./qHt(idx_ND) - p.debt_m).^(4*QuarterLength) - (1+p.rH)^(4*QuarterLength);
    stats.Spread_mean = mean(spread); % annualized
    stats.Spread_std = std(spread); % annualized
    stats.DebtServiceToOutput_mean = mean((p.debt_m + (1-p.debt_m)*p.debt_z)*bt_sum(idx_ND)./yt(idx_ND));
    stats.corr_spread_y = nancorr(spread,yt(idx_ND));
    stats.Spread_max = nanmax(spread);
    stats.RecoveryFraction_mean = nanmean(recovery_frac);
    stats.defprob1_mean = mean(defprob1t(idx_ND));
    stats.defprob1_std = std(defprob1t(idx_ND));
    stats.ConsumptionToOutput_mean = mean(ct(idx_ND)./yt(idx_ND));
    stats.volc_voly = std(ct(idx_ND))/std(yt(idx_ND));
    stats.volNXtoy_voly = std(NXt(idx_ND)./yt(idx_ND))/std(yt(idx_ND));
    stats.corr_c_y = nancorr(ct(idx_ND),yt(idx_ND));
    stats.corr_NXtoy_y = nancorr(NXt(idx_ND)./yt(idx_ND),yt(idx_ND));
    stats.bidask_ND_mean = mean(bidaskt(idx_ND));
    stats.bidask_ND_std = std(bidaskt(idx_ND));
    stats.bidask_D_mean = mean(bidaskt(idx_D));
    stats.bidask_D_std = std(bidaskt(idx_D));
    stats.corr_spread_bidask = nancorr(spread,bidaskt(idx_ND));
    stats.corr_bidask_y = nancorr(yt(idx_ND),bidaskt(idx_ND));
    TurnOver_temp = [p.prob_meet_mkt_maker_ND*reshape(1 - fracHt(idx_ND),[],1);p.prob_meet_mkt_maker_D*reshape(1 - fracHt(idx_D),[],1)];
    stats.TurnOver_mean = mean(TurnOver_temp);
    stats.TurnOver_std = std(TurnOver_temp);
    stats.TurnOver_ND_mean = mean(p.prob_meet_mkt_maker_ND*(1-fracHt(idx_ND)));
    stats.TurnOver_ND_std = std(p.prob_meet_mkt_maker_ND*(1-fracHt(idx_ND)));
    stats.TurnOver_D_mean = mean(p.prob_meet_mkt_maker_D*(1-fracHt(idx_D)));
    stats.TurnOver_D_std = std(p.prob_meet_mkt_maker_D*(1-fracHt(idx_D)));
    stats.corr_TurnOver_y_ND = nancorr(p.prob_meet_mkt_maker_ND*(1-fracHt(idx_ND)),yt(idx_ND));
    stats.corr_TurnOver_spread_ND = nancorr(p.prob_meet_mkt_maker_ND*(1-fracHt(idx_ND)),spread);
    stats.FracDebtIssue_mean = nanmean(issuefract(idx_ND));
    stats.FracDebtIssue_std = nanstd(issuefract(idx_ND));
    stats.FracDebtIssue_min = min(issuefract(idx_ND));
    stats.FracDebtIssue_max = max(issuefract(idx_ND));
end

function y = nancorr(x1,x2)
    x1 = x1(:);
    x2 = x2(:);
    idx = find(~isnan(x1) & ~isnan(x2));
    if ~isempty(idx)
        y = corr(x1(idx),x2(idx));
    else
        y = NaN;
    end
end