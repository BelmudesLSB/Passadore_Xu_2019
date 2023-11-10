%% This file compiles the mex files for the model: 

mexcuda -dynamic -DUSE_DOUBLE mexSolveModelGivenParms.cu liquidity_mex_defs.cu normaldist_mex.cu tauchen_mex.cu liquidity_vfi_mex.cu
