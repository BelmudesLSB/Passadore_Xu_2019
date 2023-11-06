function idx = sub2ind_3d_fast(SZ,idx1,idx2,idx3)

    N1 = SZ(1);
    N2 = SZ(2);
    
    idx = idx1 + (idx2-1 + (idx3-1)*N2)*N1;

end