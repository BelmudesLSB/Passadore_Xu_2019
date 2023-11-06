#include "liquidity_mex_defs.h"
#include "realtype.h"
#include "mex.h"
#include "matrix.h"
#include "liquidity_output_to_mat.h"
#include <cstdlib>

// #define DEBUG_LIQUIDITY_OUTPUT_TO_MAT

template <class T> void write_scalar_to_mat_struc(mxArray *pm, const char *fieldname, T write_val)
{

	mxArray* temp_mxArray;
	double* temp_dbl;

	temp_mxArray = mxCreateDoubleMatrix(1, 1, mxREAL);
	temp_dbl = mxGetPr(temp_mxArray);
	temp_dbl[0] = (double) write_val;
	mxSetField(pm, 0, fieldname, temp_mxArray);

}

template <class T> void write_vector_to_mat_struc(mxArray* pm, const char* fieldname, T *write_val, int num_to_write)
{

	mxArray* temp_mxArray;
	double* temp_dbl;

	temp_mxArray = mxCreateDoubleMatrix(num_to_write, 1, mxREAL);
	temp_dbl = mxGetPr(temp_mxArray);
	for (int i = 0; i < num_to_write; i++) {
		temp_dbl[i] = (double) write_val[i];
	}
	mxSetField(pm, 0, fieldname, temp_mxArray);

}

void save_parms_to_mat_struc(mxArray *pOut, parms_bsl_mod &p)
{

	const char* fieldnames[NUM_PARM_FIELDS_TO_OUTPUT];

	for (int i = 0; i < NUM_PARM_FIELDS_TO_OUTPUT; i++) {
		fieldnames[i] = (char*)mxMalloc(30);
	}

	memcpy((void*)fieldnames[0], "rra", sizeof("rra"));
	memcpy((void*)fieldnames[1], "beta", sizeof("beta"));
	/*
	memcpy((void*)fieldnames[2], "clb", sizeof("clb"));
	memcpy((void*)fieldnames[3], "ggq_mDval", sizeof("ggq_mDval"));
	memcpy((void*)fieldnames[4], "ggq_mlb", sizeof("ggq_mlb"));
	memcpy((void*)fieldnames[5], "ggq_mub", sizeof("ggq_mub"));
	memcpy((void*)fieldnames[6], "ggq_mmean", sizeof("ggq_mmean"));
	memcpy((void*)fieldnames[7], "ggq_mstd", sizeof("ggq_mstd"));
	memcpy((void*)fieldnames[8], "ggq_mmass", sizeof("ggq_mmass"));
	memcpy((void*)fieldnames[9], "errtol_bisect", sizeof("errtol_bisect"));
	memcpy((void*)fieldnames[10], "maxiters_bisect", sizeof("maxiters_bisect"));
	memcpy((void*)fieldnames[11], "rH", sizeof("rH"));
	memcpy((void*)fieldnames[12], "rL", sizeof("rL"));
	memcpy((void*)fieldnames[13], "holding_cost", sizeof("holding_cost"));
	memcpy((void*)fieldnames[14], "prob_liqshock", sizeof("prob_liqshock"));
	memcpy((void*)fieldnames[15], "prob_meet_mkt_maker_ND", sizeof("prob_meet_mkt_maker_ND"));
	memcpy((void*)fieldnames[15], "prob_meet_mkt_maker_D", sizeof("prob_meet_mkt_maker_D"));
	memcpy((void*)fieldnames[16], "bargain_power_mkt_maker_ND", sizeof("bargain_power_mkt_maker_ND"));
	memcpy((void*)fieldnames[17], "bargain_power_mkt_maker_D", sizeof("bargain_power_mkt_maker_D"));
	memcpy((void*)fieldnames[18], "reenterprob", sizeof("reenterprob"));
	memcpy((void*)fieldnames[19], "debt_z", sizeof("debt_z"));
	memcpy((void*)fieldnames[20], "debt_m", sizeof("debt_m"));
	memcpy((void*)fieldnames[21], "debt_writedown_reentry", sizeof("debt_writedown_reentry"));
	memcpy((void*)fieldnames[22], "max_def_prob", sizeof("max_def_prob"));
	memcpy((void*)fieldnames[23], "recovery_bmax", sizeof("recovery_bmax"));
	memcpy((void*)fieldnames[24], "gridsize_tfp", sizeof("gridsize_tfp"));
	memcpy((void*)fieldnames[25], "gridsize_b", sizeof("gridsize_b"));
	memcpy((void*)fieldnames[26], "bmax", sizeof("bmax"));
	memcpy((void*)fieldnames[27], "rho_tfp", sizeof("rho_tfp"));
	memcpy((void*)fieldnames[28], "sigma_tfp", sizeof("sigma_tfp"));
	memcpy((void*)fieldnames[29], "d_y", sizeof("d_y"));
	memcpy((void*)fieldnames[30], "d_yy", sizeof("d_yy"));
	memcpy((void*)fieldnames[31], "d_b", sizeof("d_b"));
	memcpy((void*)fieldnames[32], "d_bb", sizeof("d_bb"));
	memcpy((void*)fieldnames[33], "d_yb", sizeof("d_yb"));
	memcpy((void*)fieldnames[34], "grid_y_nd", sizeof("grid_y_nd"));
	memcpy((void*)fieldnames[35], "grid_y_d", sizeof("grid_y_d"));
	memcpy((void*)fieldnames[36], "grid_b", sizeof("grid_b"));
	memcpy((void*)fieldnames[37], "prob_y", sizeof("prob_y"));
	memcpy((void*)fieldnames[38], "grid_uD_yb", sizeof("grid_uD_yb"));
	memcpy((void*)fieldnames[39], "maxiters", sizeof("maxiters"));
	memcpy((void*)fieldnames[40], "errtolV", sizeof("errtolV"));
	memcpy((void*)fieldnames[41], "errtolQ", sizeof("errtolQ"));
	memcpy((void*)fieldnames[42], "NumUpdateRegimes", sizeof("NumUpdateRegimes"));
	memcpy((void*)fieldnames[43], "UpdateRegimeIter", sizeof("UpdateRegimeIter"));
	memcpy((void*)fieldnames[44], "UpdateWeightOldV", sizeof("UpdateWeightOldV"));
	memcpy((void*)fieldnames[45], "UpdateWeightOldQ", sizeof("UpdateWeightOldQ"));
	memcpy((void*)fieldnames[46], "UpdateWeightOldDefProb", sizeof("UpdateWeightOldDefProb"));
	*/

	/*
	memcpy((void*)fieldnames[47], "iters", sizeof("iters"));
	memcpy((void*)fieldnames[48], "err_CVD", sizeof("err_CVD"));
	memcpy((void*)fieldnames[49], "err_CVND", sizeof("err_CVND"));
	memcpy((void*)fieldnames[50], "err_qH_ND", sizeof("err_qH_ND"));
	memcpy((void*)fieldnames[51], "err_qL_ND", sizeof("err_qL_ND"));
	memcpy((void*)fieldnames[52], "err_qH_D", sizeof("err_qH_D"));
	memcpy((void*)fieldnames[53], "err_qL_D", sizeof("err_qL_D"));
	*/

	pOut = mxCreateStructMatrix(1,1, NUM_PARM_FIELDS_TO_OUTPUT,fieldnames);

	for (int i = 0; i < NUM_PARM_FIELDS_TO_OUTPUT; i++) {
		mxFree((void*)fieldnames[i]);
	}

	write_scalar_to_mat_struc(pOut, "rra", p.rra);
	write_scalar_to_mat_struc(pOut, "beta", p.beta);
	/*
	write_scalar_to_mat_struc(pOut, "clb", p.clb);
	write_scalar_to_mat_struc(pOut, "ggq_mDval", p.ggq_mDval);
	write_scalar_to_mat_struc(pOut, "ggq_mlb", p.ggq_mlb);
	write_scalar_to_mat_struc(pOut, "ggq_mub", p.ggq_mub);
	write_scalar_to_mat_struc(pOut, "ggq_mmean", p.ggq_mmean);
	write_scalar_to_mat_struc(pOut, "ggq_mstd", p.ggq_mstd);
	write_scalar_to_mat_struc(pOut, "ggq_mmass", p.ggq_mmass);
	write_scalar_to_mat_struc(pOut, "errtol_bisect", p.errtol_bisect);
	write_scalar_to_mat_struc(pOut, "maxiters_bisect", p.maxiters_bisect);
	write_scalar_to_mat_struc(pOut, "rH", p.rH);
	write_scalar_to_mat_struc(pOut, "rL", p.rL);
	write_scalar_to_mat_struc(pOut, "holding_cost", p.holding_cost);
	write_scalar_to_mat_struc(pOut, "prob_liqshock", p.prob_liqshock);
	write_scalar_to_mat_struc(pOut, "prob_meet_mkt_maker_ND", p.prob_meet_mkt_maker_ND);
	write_scalar_to_mat_struc(pOut, "prob_meet_mkt_maker_D", p.prob_meet_mkt_maker_D);
	write_scalar_to_mat_struc(pOut, "bargain_power_mkt_maker_ND", p.bargain_power_mkt_maker_ND);
	write_scalar_to_mat_struc(pOut, "bargain_power_mkt_maker_D", p.bargain_power_mkt_maker_D);
	write_scalar_to_mat_struc(pOut, "reenterprob", p.reenterprob);
	write_scalar_to_mat_struc(pOut, "debt_z", p.debt_z);
	write_scalar_to_mat_struc(pOut, "debt_m", p.debt_m);
	write_scalar_to_mat_struc(pOut, "debt_writedown_reentry", p.debt_writedown_reentry);
	write_scalar_to_mat_struc(pOut, "max_def_prob", p.max_def_prob);
	write_scalar_to_mat_struc(pOut, "recovery_bmax", p.recovery_bmax);
	write_scalar_to_mat_struc(pOut, "gridsize_tfp", p.gridsize_tfp);
	write_scalar_to_mat_struc(pOut, "gridsize_b", p.gridsize_b);
	write_scalar_to_mat_struc(pOut, "bmax", p.bmax);
	write_scalar_to_mat_struc(pOut, "rho_tfp", p.rho_tfp);
	write_scalar_to_mat_struc(pOut, "sigma_tfp", p.sigma_tfp);
	write_scalar_to_mat_struc(pOut, "d_y", p.d_y);
	write_scalar_to_mat_struc(pOut, "d_yy", p.d_yy);
	write_scalar_to_mat_struc(pOut, "d_b", p.d_b);
	write_scalar_to_mat_struc(pOut, "d_bb", p.d_bb);
	write_scalar_to_mat_struc(pOut, "d_yb", p.d_yb);
	write_vector_to_mat_struc(pOut, "grid_y_nd", p.grid_y_nd, p.gridsize_tfp);
	write_vector_to_mat_struc(pOut, "grid_y_d", p.grid_y_d, p.gridsize_tfp);
	write_vector_to_mat_struc(pOut, "grid_b", p.grid_b, p.gridsize_b);
	write_vector_to_mat_struc(pOut, "prob_y", p.prob_y, p.gridsize_tfp* p.gridsize_tfp);
	write_vector_to_mat_struc(pOut, "grid_uD_yb", p.grid_uD_yb, p.gridsize_tfp * p.gridsize_b);
	write_scalar_to_mat_struc(pOut, "maxiters", p.maxiters);
	write_scalar_to_mat_struc(pOut, "errtolV", p.errtolV);
	write_scalar_to_mat_struc(pOut, "errtolQ", p.errtolQ);
	write_scalar_to_mat_struc(pOut, "NumUpdateRegimes", p.NumUpdateRegimes);
	write_vector_to_mat_struc(pOut, "UpdateRegimeIter", p.UpdateRegimeIter, p.NumUpdateRegimes - 1);
	write_vector_to_mat_struc(pOut, "UpdateWeightOldV", p.UpdateWeightOldV, p.NumUpdateRegimes);
	write_vector_to_mat_struc(pOut, "UpdateWeightOldQ", p.UpdateWeightOldQ, p.NumUpdateRegimes);
	write_vector_to_mat_struc(pOut, "UpdateWeightOldDefProb", p.UpdateWeightOldDefProb, p.NumUpdateRegimes);
	*/

}