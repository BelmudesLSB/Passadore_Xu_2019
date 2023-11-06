function SPARSE_TRANSITION_PROB_YBD = ...
    GenSparseTransitionMatrix(p, defprob, idx_bnext, proby_cutoff)
% this algorithm uses two approximations to speed things up:
% 1. cutoff for the transition f(y'|y)
% 2. switch off the randomization shocks eps' when it comes to policies

    if proby_cutoff<0 || proby_cutoff>1; error('proby_cutoff should be between 0 and 1'); end
    if any(size(defprob)~=[p.gridsize_tfp p.gridsize_b]); error('defprob does not have the correct dimensions'); end
    if any(size(idx_bnext)~=[p.gridsize_tfp p.gridsize_b]); error('idx_bnext does not have the correct dimensions'); end
    
    
    % gridsize_yb = p.gridsize_tfp*p.gridsize_b;
    gridsize_ybd = p.gridsize_tfp*p.gridsize_b*2;
    SZ_YBD = [p.gridsize_tfp p.gridsize_b 2];
    
    %% 1. precompute
    
    % sparse transition matrix for y
    num_ynext = zeros(p.gridsize_tfp,1);
    idx_ynext = cell(p.gridsize_tfp,1);
    prob_ynext = cell(p.gridsize_tfp,1);
    prob_ynext_not_reenter = cell(p.gridsize_tfp,1);
    prob_ynext_reenter = cell(p.gridsize_tfp,1);
    for idx_y = 1:p.gridsize_tfp
        idx_ynext{idx_y} = find(p.prob_y(idx_y,:)>proby_cutoff);
        num_ynext(idx_y) = numel(idx_ynext{idx_y});
        prob_ynext{idx_y} = reshape(p.prob_y(idx_y,idx_ynext{idx_y})/sum(p.prob_y(idx_y,idx_ynext{idx_y})),[],1);
        prob_ynext_not_reenter{idx_y} = (1 - p.reenterprob)*prob_ynext{idx_y};
        prob_ynext_reenter{idx_y} = p.reenterprob*prob_ynext{idx_y};
    end
    
    % reenter
    idx_b_reenter = zeros(p.gridsize_b,1);
    for ib = 1:p.gridsize_b
        [~,idx_b_reenter(ib)] = min(abs(p.grid_b - min((1-p.debt_writedown_reentry)*p.grid_b(ib),p.recovery_bmax)));
    end
    
    % indices
    SPARSE_I = zeros(5*sum(num_ynext),1); % from
    SPARSE_J = zeros(5*sum(num_ynext),1); % to
    SPARSE_V = zeros(5*sum(num_ynext),1); % prob(from->to)
    
    %% 2. transition matrix, in default (DEF=1)
    
    IDX_SPARSE_END = 0;
    
    for idx_b = 1:p.gridsize_b
        for idx_y = 1:p.gridsize_tfp
            % not reenter
            IDX_SPARSE_START = IDX_SPARSE_END + 1;
            IDX_SPARSE_END = IDX_SPARSE_END + num_ynext(idx_y);
            SPARSE_I(IDX_SPARSE_START:IDX_SPARSE_END) = sub2ind_3d_fast(SZ_YBD, idx_y, idx_b, 1);
            SPARSE_J(IDX_SPARSE_START:IDX_SPARSE_END) = sub2ind_3d_fast(SZ_YBD, idx_ynext{idx_y}, idx_b, 1);
            SPARSE_V(IDX_SPARSE_START:IDX_SPARSE_END) = prob_ynext_not_reenter{idx_y};
            % reenter, not default
            IDX_SPARSE_START = IDX_SPARSE_END + 1;
            IDX_SPARSE_END = IDX_SPARSE_END + num_ynext(idx_y);
            SPARSE_I(IDX_SPARSE_START:IDX_SPARSE_END) = sub2ind_3d_fast(SZ_YBD, idx_y, idx_b, 1);
            SPARSE_J(IDX_SPARSE_START:IDX_SPARSE_END) = sub2ind_3d_fast(SZ_YBD, idx_ynext{idx_y}, idx_b_reenter(idx_b), 2);
            SPARSE_V(IDX_SPARSE_START:IDX_SPARSE_END) = prob_ynext_reenter{idx_y}.*(1-defprob(idx_ynext{idx_y},idx_b_reenter(idx_b)));
            % reenter, default
            IDX_SPARSE_START = IDX_SPARSE_END + 1;
            IDX_SPARSE_END = IDX_SPARSE_END + num_ynext(idx_y);
            SPARSE_I(IDX_SPARSE_START:IDX_SPARSE_END) = sub2ind_3d_fast(SZ_YBD, idx_y, idx_b, 1);
            SPARSE_J(IDX_SPARSE_START:IDX_SPARSE_END) = sub2ind_3d_fast(SZ_YBD, idx_ynext{idx_y}, idx_b_reenter(idx_b), 1);
            SPARSE_V(IDX_SPARSE_START:IDX_SPARSE_END) = prob_ynext_reenter{idx_y}.*defprob(idx_ynext{idx_y},idx_b_reenter(idx_b));
        end
    end
    
    
    %% 3. transition matrix, not in default (DEF=2)
    
    for idx_b = 1:p.gridsize_b
        for idx_y = 1:p.gridsize_tfp
            % default
            IDX_SPARSE_START = IDX_SPARSE_END + 1;
            IDX_SPARSE_END = IDX_SPARSE_END + num_ynext(idx_y);
            SPARSE_I(IDX_SPARSE_START:IDX_SPARSE_END) = sub2ind_3d_fast(SZ_YBD, idx_y, idx_b, 2);
            SPARSE_J(IDX_SPARSE_START:IDX_SPARSE_END) = sub2ind_3d_fast(SZ_YBD, idx_ynext{idx_y}, idx_bnext(idx_y,idx_b), 1);
            SPARSE_V(IDX_SPARSE_START:IDX_SPARSE_END) = prob_ynext{idx_y}.*defprob(idx_ynext{idx_y},idx_bnext(idx_y,idx_b));
            % not default
            IDX_SPARSE_START = IDX_SPARSE_END + 1;
            IDX_SPARSE_END = IDX_SPARSE_END + num_ynext(idx_y);
            SPARSE_I(IDX_SPARSE_START:IDX_SPARSE_END) = sub2ind_3d_fast(SZ_YBD, idx_y, idx_b, 2);
            SPARSE_J(IDX_SPARSE_START:IDX_SPARSE_END) = sub2ind_3d_fast(SZ_YBD, idx_ynext{idx_y}, idx_bnext(idx_y,idx_b), 2);
            SPARSE_V(IDX_SPARSE_START:IDX_SPARSE_END) = prob_ynext{idx_y}.*(1-defprob(idx_ynext{idx_y},idx_bnext(idx_y,idx_b)));
        end
    end
    
    %% 4. create transition matrix
    
    SPARSE_TRANSITION_PROB_YBD = sparse(SPARSE_I,SPARSE_J,SPARSE_V,gridsize_ybd,gridsize_ybd);

end