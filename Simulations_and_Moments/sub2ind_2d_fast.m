function idx = sub2ind_2d_fast(SZ,idx1,idx2)

    N1 = SZ(1);
    % N2 = SZ(2);
    
    idx = idx1 + N1*(idx2 - 1);

end