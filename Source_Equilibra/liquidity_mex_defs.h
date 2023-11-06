#ifndef LIQUIDITY_MEX_DEFS_H
#define LIQUIDITY_MEX_DEFS_H

#include "realtype.h"
#include "mex.h"

#ifndef MAX_GRIDSIZE_Y
#define MAX_GRIDSIZE_Y 256
#endif

#ifndef MAX_GRIDSIZE_B
#define MAX_GRIDSIZE_B 512
#endif

REAL_TYPE CalcAutarkyOutput(REAL_TYPE yND, REAL_TYPE bdef, REAL_TYPE d_y, REAL_TYPE d_yy, REAL_TYPE d_b, REAL_TYPE d_bb, REAL_TYPE d_yb);

int ReadScalarFromMatlabStruc(double& ReadVal, const mxArray* mxPtr, const char* fieldname);
int ReadVectorFromMatlabStruc(double* ReadVal, const mxArray* mxPtr, int Nread, const char* fieldname);

class parms_bsl_mod
{
public:
	// preferences
	REAL_TYPE rra, beta;

	// numerical parameters
	REAL_TYPE clb, ggq_mDval, ggq_mlb, ggq_mub, ggq_mmean, ggq_mstd, ggq_mmass, errtol_bisect;
	int maxiters_bisect;

	// liquidity discounts
	REAL_TYPE rH, rL, holding_cost;

	// liquidity transitions
	REAL_TYPE prob_liqshock;
	REAL_TYPE prob_meet_mkt_maker_ND, prob_meet_mkt_maker_D;

	// bargaining powers
	REAL_TYPE bargain_power_mkt_maker_ND, bargain_power_mkt_maker_D;

	// debt market
	REAL_TYPE reenterprob, debt_z, debt_m;
	REAL_TYPE debt_writedown_reentry;
	REAL_TYPE recovery_bmax;
	REAL_TYPE max_def_prob;

	// grid dimensions
	REAL_TYPE tfp_lb, tfp_ub;
	REAL_TYPE conv_chk_tfp_lb, conv_chk_tfp_ub, conv_chk_b_ub;
	int gridsize_tfp, gridsize_b;
	REAL_TYPE bmax;

	// output
	REAL_TYPE rho_tfp, sigma_tfp;
	REAL_TYPE d_y, d_yy, d_b, d_bb, d_yb; // default costs

	// grids
	REAL_TYPE *grid_y_nd, *grid_y_d, *grid_b;
	REAL_TYPE *prob_y;
	REAL_TYPE *grid_uD_yb; 

	// convergence
	int maxiters;
	REAL_TYPE errtolV, errtolQ;
	int NumUpdateRegimes; // e.g. 3
	int *UpdateRegimeIter; // e.g. 100/500
	REAL_TYPE *UpdateWeightOldV, *UpdateWeightOldQ, *UpdateWeightOldDefProb; // 0/0.5/0.9

	int checkparms();

	int readparms(const mxArray* mxPtr); // defined in liquidity_io_matlab.cpp

	void matlab_display_cpu_parms();

	void initDeviceConstants(); // defined in vfi

	void clean();

};

#endif