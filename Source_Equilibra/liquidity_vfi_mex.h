#ifndef LIQUIDITY_VFI_MEX_H
#define LIQUIDITY_VFI_MEX_H

#include "realtype.h"
#include "liquidity_mex_defs.h"

__device__ void ggq_topdown(REAL_TYPE& CENDupdate, REAL_TYPE& qHupdate, REAL_TYPE& qLupdate, REAL_TYPE& def_prob_update, REAL_TYPE& mdeftresh,
	int idx_y, int idx_b, REAL_TYPE& VD, REAL_TYPE* d_qH_ND, REAL_TYPE* d_qL_ND, REAL_TYPE* d_qH_D, REAL_TYPE* d_qL_D, REAL_TYPE* d_CVND, REAL_TYPE* d_def_prob);

__device__ REAL_TYPE bisect_zero(REAL_TYPE a, REAL_TYPE b, REAL_TYPE c1, REAL_TYPE c2, REAL_TYPE W1_minus_W2);
__device__ REAL_TYPE m_root_fun(REAL_TYPE m, REAL_TYPE U1, REAL_TYPE U2, REAL_TYPE W1_minus_W2);
__device__ REAL_TYPE CENDupdatefun(REAL_TYPE m, REAL_TYPE U, REAL_TYPE W);
__device__ REAL_TYPE gauss_legendre_CENDupdate(REAL_TYPE m1, REAL_TYPE m2, REAL_TYPE U, REAL_TYPE W);

__global__ void vfi_iterate_policy(int* d_idx_bchoice, REAL_TYPE* d_qH_ND, REAL_TYPE* d_def_prob, REAL_TYPE* d_CVND);

__global__ void vfi_iterate(REAL_TYPE *d_VD, REAL_TYPE *d_VNDupdate, REAL_TYPE *d_qH_ND_update, REAL_TYPE *d_qL_ND_update, 
	REAL_TYPE *d_defprob_update, REAL_TYPE *d_defthresh,
	REAL_TYPE *d_CVD, REAL_TYPE *d_qH_ND, REAL_TYPE *d_qL_ND, 
	REAL_TYPE *d_qH_D, REAL_TYPE *d_qL_D, REAL_TYPE *d_CVND, REAL_TYPE *d_uD_yb, REAL_TYPE *d_defprob);

__global__ void vfi_interpolate(REAL_TYPE *ceCi, REAL_TYPE *qH_NDi, REAL_TYPE *qL_NDi,
	REAL_TYPE *ceC, REAL_TYPE *qH_ND, REAL_TYPE *qL_ND);

__global__ void vfi_update1(REAL_TYPE *d_prob_y, REAL_TYPE *d_EVD, REAL_TYPE *d_ceC, 
	REAL_TYPE *d_qH_ND, REAL_TYPE *d_qL_ND, REAL_TYPE *d_qH_D1, REAL_TYPE *d_qL_D1, REAL_TYPE *d_Edefprob,
	REAL_TYPE *d_VD, REAL_TYPE *d_VNDupdate, REAL_TYPE *d_qH_NDupdate, REAL_TYPE *d_qL_NDupdate, REAL_TYPE *d_defprob_update, 
	REAL_TYPE *d_qH_D0, REAL_TYPE *d_qL_D0);

__global__ void vfi_update2(REAL_TYPE *d_CVD, REAL_TYPE *d_CVDhaircut, REAL_TYPE *d_CVNDhaircut,
	REAL_TYPE *d_qH_D, REAL_TYPE *d_qL_D, REAL_TYPE *d_qH_Dhaircut, REAL_TYPE *d_qL_Dhaircut,
	REAL_TYPE *d_qH_NDhaircut, REAL_TYPE *d_qL_NDhaircut);

__global__ void update_compute_errors(REAL_TYPE *d_x, REAL_TYPE *d_y, REAL_TYPE wx);

void vfi(parms_bsl_mod &p, REAL_TYPE wOldV, REAL_TYPE wOldQ, REAL_TYPE wOldDefProb, REAL_TYPE *d_prob_y, REAL_TYPE *d_uD_yb, 
	REAL_TYPE &err_CVD, REAL_TYPE &err_CVND, REAL_TYPE &err_qH_ND, REAL_TYPE &err_qL_ND, REAL_TYPE &err_qH_D, REAL_TYPE &err_qL_D, REAL_TYPE &err_defprob,
	REAL_TYPE *d_CVD, REAL_TYPE *d_CVND, REAL_TYPE *d_qH_ND, REAL_TYPE *d_qL_ND, REAL_TYPE *d_qH_D, REAL_TYPE *d_qL_D, REAL_TYPE *d_defprob, REAL_TYPE *d_VD, 
	REAL_TYPE *d_VNDupdate, REAL_TYPE *d_qH_ND_update, REAL_TYPE *d_qL_ND_update, REAL_TYPE *d_defprob_update, REAL_TYPE *d_defthresh,
	REAL_TYPE *d_EVD, REAL_TYPE *d_CVNDnew,
	REAL_TYPE *d_qH_NDnew, REAL_TYPE *d_qL_NDnew, REAL_TYPE *d_EqH_D, REAL_TYPE *d_EqL_D, REAL_TYPE *d_defprob_new,
	REAL_TYPE *d_CVNDi, REAL_TYPE *d_qH_NDi, REAL_TYPE *d_qL_NDi,
	REAL_TYPE *d_CVDnew, REAL_TYPE *d_qH_Dnew, REAL_TYPE *d_qL_Dnew);

int SolveModel(REAL_TYPE* h_CVD, REAL_TYPE* h_CVND, REAL_TYPE* h_qH_ND, REAL_TYPE* h_qL_ND, REAL_TYPE* h_qH_D, REAL_TYPE* h_qL_D,
	REAL_TYPE* h_defprob, REAL_TYPE* h_defthresh, int* h_idx_bchoice,
	int& iter, REAL_TYPE& err_qH_ND, REAL_TYPE& err_qL_ND, REAL_TYPE& err_qH_D, REAL_TYPE& err_qL_D, REAL_TYPE& err_defprob,
	REAL_TYPE& err_CVD, REAL_TYPE& err_CVND,
	parms_bsl_mod& p, int useDevice = 0, bool bDisplayProgress = true, int DisplayInterval = 25);

__global__ void initValueFuns(REAL_TYPE *d_CVD, REAL_TYPE *d_CVND, REAL_TYPE *d_qH_ND, REAL_TYPE *d_qL_ND, REAL_TYPE *d_qH_D, REAL_TYPE *d_qL_D,
	REAL_TYPE *d_defprob);

#endif