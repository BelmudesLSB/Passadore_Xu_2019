#include "liquidity_mex_defs.h"
#include "realtype.h"
#include "tauchen.h"
#include <cstdlib>
#include <cmath>
#include "normaldist_mex.h"
#include "mex.h"

// #define DEBUG_IDENTITY_PROB

/*
// preferences
__constant__ REAL_TYPE RRA;
__constant__ REAL_TYPE BETA;
__constant__ REAL_TYPE ONE_MINUS_RRA; // 1-RRA
__constant__ REAL_TYPE U_SCALING; // 1/(1-RRA)

// numerical parameters
__constant__ REAL_TYPE C_LB;
__constant__ REAL_TYPE GGQ_MLB;
__constant__ REAL_TYPE GGQ_MUB;
__constant__ REAL_TYPE GGQ_MMEAN;
__constant__ REAL_TYPE GGQ_MSTD;
__constant__ REAL_TYPE GGQ_MMASS;
__constant__ REAL_TYPE GGQ_MDEF;
__constant__ REAL_TYPE ERRTOL_BISECT;
__constant__ REAL_TYPE MAX_DEF_PROB;
__constant__ int MAXITERS_BISECT;
__constant__ int CONV_CHK_TFP_LB_IDX;
__constant__ int CONV_CHK_TFP_UB_IDX;

// liquidity discounts
__constant__ REAL_TYPE RH;
__constant__ REAL_TYPE RL;
__constant__ REAL_TYPE HOLDING_COST;

//// effective transitions

__constant__ REAL_TYPE PROB_HL;
__constant__ REAL_TYPE PROB_LH_ND;
__constant__ REAL_TYPE PROB_LH_D;

// Note:
// 1. H->L is Prob(get liquidity shock)
// 2. L->H is Prob(meet mkt maker)*(1-(mkt maker barg power))

// debt market
__constant__ REAL_TYPE REENTERPROB;
__constant__ REAL_TYPE DEBT_Z;
__constant__ REAL_TYPE DEBT_M;

// grid dimensions
__constant__ int GRIDSIZE_Y;
__constant__ int GRIDSIZE_B;

// grids
__constant__ REAL_TYPE GRID_Y_ND[MAX_GRIDSIZE_Y];
__constant__ REAL_TYPE GRID_B[MAX_GRIDSIZE_B];

// interpolation
__constant__ int GRID_B_REENTER_IDX[MAX_GRIDSIZE_B];
__constant__ REAL_TYPE GRID_RECOVERY_FRAC[MAX_GRIDSIZE_B];

*/

REAL_TYPE CalcAutarkyOutput(REAL_TYPE yND, REAL_TYPE bdef, REAL_TYPE d_y, REAL_TYPE d_yy, REAL_TYPE d_b, REAL_TYPE d_bb, REAL_TYPE d_yb)
{
	REAL_TYPE yD = yND - min(yND, max(0.0, d_y*yND + d_yy*yND*yND + d_b*bdef + d_bb*bdef*bdef + d_yb*yND*bdef));
	return yD;
}

void parms_bsl_mod::clean()
{
	if (grid_y_nd)
	{
		delete[] grid_y_nd;
		grid_y_nd = NULL;
	}
	if (grid_y_d)
	{
		delete[] grid_y_d;
		grid_y_d = NULL;
	}
	if (grid_b)
	{
		delete[] grid_b;
		grid_b = NULL;
	}
	if (prob_y)
	{
		delete[] prob_y;
		prob_y = NULL;
	}
	if (grid_uD_yb)
	{
		delete[] grid_uD_yb;
		grid_uD_yb = NULL;
	}
	if (UpdateRegimeIter)
	{
		delete[] UpdateRegimeIter;
		UpdateRegimeIter = NULL;
	}
	if (UpdateWeightOldV)
	{
		delete[] UpdateWeightOldV;
		UpdateWeightOldV = NULL;
	}
	if (UpdateWeightOldQ)
	{
		delete[] UpdateWeightOldQ;
		UpdateWeightOldQ = NULL;
	}
	if (UpdateWeightOldDefProb)
	{
		delete[] UpdateWeightOldDefProb;
		UpdateWeightOldDefProb = NULL;
	}
}

int parms_bsl_mod::checkparms()
{
	bool bParmsAllOK = true;
	if (rra <= 0)
	{
		mexPrintf("warning: rra<=0.\n");
		bParmsAllOK = false;
	}
	if (beta<0 || beta>=1)
	{
		mexPrintf("warning: beta<0 || beta>=1.\n");
		bParmsAllOK = false;
	}
	if (clb <= 0)
	{
		mexPrintf("warning: clb<=0.\n");
		bParmsAllOK = false;
	}
	if (ggq_mlb >= ggq_mub)
	{
		mexPrintf("warning: ggq_mlb >= ggq_mub.\n");
		bParmsAllOK = false;
	}
	if (ggq_mDval < ggq_mlb || ggq_mDval > ggq_mub)
	{
		mexPrintf("warning: ggq_mDval < ggq_mlb || ggq_mDval > ggq_mub.\n");
		bParmsAllOK = false;
	}
	if (ggq_mmean < ggq_mlb || ggq_mmean > ggq_mub)
	{
		mexPrintf("warning: ggq_mmean < ggq_mlb || ggq_mmean > ggq_mub.\n");
		bParmsAllOK = false;
	}
	if (ggq_mstd <= 0 )
	{
		mexPrintf("warning: ggq_mstd <= 0.\n");
		bParmsAllOK = false;
	}
	if (errtol_bisect <= 0)
	{
		mexPrintf("warning: errtol_bisect <= 0.\n");
		bParmsAllOK = false;
	}
	if (maxiters_bisect <= 0)
	{
		mexPrintf("warning: maxiters_bisect <= 0.\n");
		bParmsAllOK = false;
	}
	if (rH > rL)
	{
		mexPrintf("warning: rH > rL.\n");
		bParmsAllOK = false;
	}
	if (rH < 0 || rL < 0)
	{
		mexPrintf("warning: rH < 0 || rL < 0.\n");
		bParmsAllOK = false;
	}
	if (prob_liqshock<0 || prob_liqshock>1 ||
		prob_meet_mkt_maker_ND<0 || prob_meet_mkt_maker_ND>1 ||
		prob_meet_mkt_maker_D < 0 || prob_meet_mkt_maker_D>1 ||
		bargain_power_mkt_maker_ND<0 || bargain_power_mkt_maker_ND>1 ||
		bargain_power_mkt_maker_D<0 || bargain_power_mkt_maker_D>1)
	{
		mexPrintf("warning: probabilities from liquidity transitions not valid.\n");
		bParmsAllOK = false;
	}

	if (reenterprob < 0 || reenterprob>1)
	{
		mexPrintf("warning: reenterprob not valid.\n");
		bParmsAllOK = false;
	}

	if (debt_z < 0 || debt_m<0 || debt_m>1 ||  
		debt_writedown_reentry<0 || debt_writedown_reentry>1 || max_def_prob<0 || max_def_prob>1 || recovery_bmax<0)
	{
		mexPrintf("warning: invalid debt parameters.\n");
		bParmsAllOK = false;
	}

	if (gridsize_tfp<2 || gridsize_b<2)
	{
		mexPrintf("warning: too few grid points.\n");
		bParmsAllOK = false;
	}

	if (gridsize_tfp > MAX_GRIDSIZE_Y || gridsize_b>MAX_GRIDSIZE_B)
	{
		mexPrintf("warning: grid size exceeds maximum allowed (%d for y, %d for b).\n", MAX_GRIDSIZE_Y, MAX_GRIDSIZE_B);
		bParmsAllOK = false;
	}

	if (bmax<=0)
	{
		mexPrintf("warning: bmax<=0.\n");
		bParmsAllOK = false;
	}

	if (rho_tfp<=0 || rho_tfp>=1 || sigma_tfp<=0 || tfp_lb>=0 || tfp_ub<=0 || conv_chk_tfp_lb>=conv_chk_tfp_ub)
	{
		mexPrintf("warning: invalid tfp dynamics.\n");
		bParmsAllOK = false;
	}

	bool bValidUpdateWeights = true;
	for (int i = 0; i < NumUpdateRegimes; i++)
	{
		if (UpdateWeightOldV[i]<0 || UpdateWeightOldV[i]>1 ||
			UpdateWeightOldQ[i]<0 || UpdateWeightOldQ[i]>1 ||
			UpdateWeightOldDefProb[i]<0 || UpdateWeightOldDefProb[i]>1)
		{
			bValidUpdateWeights = false;
			break;
		}
	}
	if (!bValidUpdateWeights)
	{
		mexPrintf("warning: invalid update weights detected.\n");
		bParmsAllOK = false;
	}
	
	if (bParmsAllOK)
	{
		mexPrintf("Parameters satisfy checks.\n");
		return EXIT_SUCCESS;
	}
	else
	{
		mexPrintf("Parameters do not satisfy checks.\n");
		return EXIT_FAILURE;
	}

}

int ReadScalarFromMatlabStruc(double& ReadVal, const mxArray* mxPtr, const char* fieldname) {

	if (mxGetFieldNumber(mxPtr, fieldname) == -1)
	{
		// field not found
		return(EXIT_FAILURE);
	}
	else
	{
		mxArray* ptr_mxarr = mxGetField(mxPtr, 0, fieldname);

		if (ptr_mxarr != NULL)
		{
			double* ptr_real = mxGetPr(ptr_mxarr);

			ReadVal = *ptr_real;

			return(EXIT_SUCCESS);
		}
		else
		{
			return(EXIT_FAILURE);
		}
	}

}

int ReadVectorFromMatlabStruc(double* ReadVal, const mxArray* mxPtr, int Nread, const char* fieldname) {

	if (mxGetFieldNumber(mxPtr, fieldname) == -1)
	{
		// field not found
		return(EXIT_FAILURE);
	}
	else
	{

		mxArray* ptr_mxarr = mxGetField(mxPtr, 0, fieldname);

		if (ptr_mxarr != NULL)
		{

			int sz = mxGetNumberOfElements(ptr_mxarr);

			if (sz == Nread)
			{
				double* ptr_real = mxGetPr(ptr_mxarr);

				for (int i = 0; i < Nread; i++)
				{
					ReadVal[i] = ptr_real[i];
				}

				return(EXIT_SUCCESS);
			}
			else {
				return(EXIT_FAILURE);
			}


		}
		else
		{
			return(EXIT_FAILURE);
		}
	}

}

int parms_bsl_mod::readparms(const mxArray* mxPtr) {

	double tempdbl;

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "rra") == EXIT_SUCCESS)
	{
		rra = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter rra\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "beta") == EXIT_SUCCESS)
	{
		beta = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter beta\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "clb") == EXIT_SUCCESS)
	{
		clb = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter clb\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "ggq_mDval") == EXIT_SUCCESS)
	{
		ggq_mDval = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter ggq_mDval\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "ggq_mlb") == EXIT_SUCCESS)
	{
		ggq_mlb = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter ggq_mlb\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "ggq_mub") == EXIT_SUCCESS)
	{
		ggq_mub = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter ggq_mub\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "ggq_mmean") == EXIT_SUCCESS)
	{
		ggq_mmean = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter ggq_mmean\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "ggq_mstd") == EXIT_SUCCESS)
	{
		ggq_mstd = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter ggq_mstd\n");
		return(EXIT_FAILURE);
	}

	ggq_mmass = normcdf(ggq_mub - ggq_mmean, ggq_mstd) - normcdf(ggq_mlb - ggq_mmean, ggq_mstd);

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "max_def_prob") == EXIT_SUCCESS)
	{
		max_def_prob = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter max_def_prob\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "d_y") == EXIT_SUCCESS)
	{
		d_y = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter d_y\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "d_yy") == EXIT_SUCCESS)
	{
		d_yy = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter d_yy\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "d_b") == EXIT_SUCCESS)
	{
		d_b = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter d_b\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "d_bb") == EXIT_SUCCESS)
	{
		d_bb = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter d_bb\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "d_yb") == EXIT_SUCCESS)
	{
		d_yb = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter d_yb\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "errtol_bisect") == EXIT_SUCCESS)
	{
		errtol_bisect = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter errtol_bisect\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "maxiters_bisect") == EXIT_SUCCESS)
	{
		maxiters_bisect = (int)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter maxiters_bisect\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "rH") == EXIT_SUCCESS)
	{
		rH = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter rH\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "rL") == EXIT_SUCCESS)
	{
		rL = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter rL\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "holding_cost") == EXIT_SUCCESS)
	{
		holding_cost = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter holding_cost\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "prob_liqshock") == EXIT_SUCCESS)
	{
		prob_liqshock = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter prob_liqshock\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "prob_meet_mkt_maker_ND") == EXIT_SUCCESS)
	{
		prob_meet_mkt_maker_ND = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter prob_meet_mkt_maker_ND\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "prob_meet_mkt_maker_D") == EXIT_SUCCESS)
	{
		prob_meet_mkt_maker_D = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter prob_meet_mkt_maker_D\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "bargain_power_mkt_maker_ND") == EXIT_SUCCESS)
	{
		bargain_power_mkt_maker_ND = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter bargain_power_mkt_maker_ND\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "bargain_power_mkt_maker_D") == EXIT_SUCCESS)
	{
		bargain_power_mkt_maker_D = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter bargain_power_mkt_maker_D\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "reenterprob") == EXIT_SUCCESS)
	{
		reenterprob = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter reenterprob\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "debt_z") == EXIT_SUCCESS)
	{
		debt_z = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter debt_z\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "debt_m") == EXIT_SUCCESS)
	{
		debt_m = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter debt_m\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "debt_writedown_reentry") == EXIT_SUCCESS)
	{
		debt_writedown_reentry = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter debt_writedown_reentry\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "recovery_bmax") == EXIT_SUCCESS)
	{
		recovery_bmax = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter recovery_bmax\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "gridsize_tfp") == EXIT_SUCCESS)
	{
		gridsize_tfp = (int)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter gridsize_tfp\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "gridsize_b") == EXIT_SUCCESS)
	{
		gridsize_b = (int)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter gridsize_b\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "bmax") == EXIT_SUCCESS)
	{
		bmax = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter bmax\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "sigma_tfp") == EXIT_SUCCESS)
	{
		sigma_tfp = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter sigma_tfp\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "rho_tfp") == EXIT_SUCCESS)
	{
		rho_tfp = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter rho_tfp\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "tfp_lb") == EXIT_SUCCESS)
	{
		tfp_lb = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter tfp_lb\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "tfp_ub") == EXIT_SUCCESS)
	{
		tfp_ub = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter tfp_ub\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "conv_chk_tfp_lb") == EXIT_SUCCESS)
	{
		conv_chk_tfp_lb = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter conv_chk_tfp_lb\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "conv_chk_tfp_ub") == EXIT_SUCCESS)
	{
		conv_chk_tfp_ub = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter conv_chk_tfp_ub\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "conv_chk_b_ub") == EXIT_SUCCESS)
	{
		conv_chk_b_ub = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter conv_chk_b_ub\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "maxiters") == EXIT_SUCCESS)
	{
		maxiters = (int)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter maxiters\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "errtolV") == EXIT_SUCCESS)
	{
		errtolV = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter errtolV\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "errtolQ") == EXIT_SUCCESS)
	{
		errtolQ = (REAL_TYPE)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter errtolQ\n");
		return(EXIT_FAILURE);
	}

	if (ReadScalarFromMatlabStruc(tempdbl, mxPtr, "NumUpdateRegimes") == EXIT_SUCCESS)
	{
		NumUpdateRegimes = (int)tempdbl;
	}
	else
	{
		mexPrintf("Unable to read parameter NumUpdateRegimes\n");
		return(EXIT_FAILURE);
	}

	double* tempdbl_vec;

	tempdbl_vec = new double[NumUpdateRegimes - 1];
	if (ReadVectorFromMatlabStruc(tempdbl_vec, mxPtr, NumUpdateRegimes - 1, "UpdateRegimeIter") == EXIT_SUCCESS)
	{
		UpdateRegimeIter = new int[NumUpdateRegimes - 1];
		for (int i = 0; i < NumUpdateRegimes - 1; i++)
		{
			UpdateRegimeIter[i] = (int)tempdbl_vec[i];
		}
		delete[] tempdbl_vec;
	}
	else
	{
		delete[] tempdbl_vec;
		mexPrintf("Unable to read parameter UpdateRegimeIter\n");
		return(EXIT_FAILURE);
	}

	tempdbl_vec = new double[NumUpdateRegimes];
	if (ReadVectorFromMatlabStruc(tempdbl_vec, mxPtr, NumUpdateRegimes, "UpdateWeightOldV") == EXIT_SUCCESS)
	{
		UpdateWeightOldV = new REAL_TYPE[NumUpdateRegimes];
		for (int i = 0; i < NumUpdateRegimes; i++)
		{
			UpdateWeightOldV[i] = (REAL_TYPE)tempdbl_vec[i];
		}
		delete[] tempdbl_vec;
	}
	else
	{
		delete[] tempdbl_vec;
		mexPrintf("Unable to read parameter UpdateWeightOldV\n");
		return(EXIT_FAILURE);
	}

	tempdbl_vec = new double[NumUpdateRegimes];
	if (ReadVectorFromMatlabStruc(tempdbl_vec, mxPtr, NumUpdateRegimes, "UpdateWeightOldQ") == EXIT_SUCCESS)
	{
		UpdateWeightOldQ = new REAL_TYPE[NumUpdateRegimes];
		for (int i = 0; i < NumUpdateRegimes; i++)
		{
			UpdateWeightOldQ[i] = (REAL_TYPE)tempdbl_vec[i];
		}
		delete[] tempdbl_vec;
	}
	else
	{
		delete[] tempdbl_vec;
		mexPrintf("Unable to read parameter UpdateWeightOldQ\n");
		return(EXIT_FAILURE);
	}

	tempdbl_vec = new double[NumUpdateRegimes];
	if (ReadVectorFromMatlabStruc(tempdbl_vec, mxPtr, NumUpdateRegimes, "UpdateWeightOldDefProb") == EXIT_SUCCESS)
	{
		UpdateWeightOldDefProb = new REAL_TYPE[NumUpdateRegimes];
		for (int i = 0; i < NumUpdateRegimes; i++)
		{
			UpdateWeightOldDefProb[i] = (REAL_TYPE)tempdbl_vec[i];
		}
		delete[] tempdbl_vec;
	}
	else
	{
		delete[] tempdbl_vec;
		mexPrintf("Unable to read parameter UpdateWeightOldDefProb\n");
		return(EXIT_FAILURE);
	}

	// set up the grid for b
	grid_b = new REAL_TYPE[gridsize_b];
	REAL_TYPE db = bmax / (REAL_TYPE)(gridsize_b - 1);
	for (int i = 0; i < gridsize_b; i++)
	{
		grid_b[i] = db * (REAL_TYPE)i;
	}

	// use tauchen
	REAL_TYPE* z, * prob_z;
	z = new REAL_TYPE[gridsize_tfp];
	prob_z = new REAL_TYPE[gridsize_tfp * gridsize_tfp];

	if (tauchen(z, prob_z, rho_tfp, sigma_tfp, gridsize_tfp, tfp_lb, tfp_ub) == EXIT_SUCCESS)
	{
		grid_y_nd = new REAL_TYPE[gridsize_tfp];
		grid_y_d = new REAL_TYPE[gridsize_tfp * gridsize_b];
		grid_uD_yb = new REAL_TYPE[gridsize_tfp * gridsize_b];
		prob_y = new REAL_TYPE[gridsize_tfp * gridsize_tfp];
		for (int i = 0; i < gridsize_tfp; i++)
		{
			grid_y_nd[i] = exp(z[i]);
			for (int j = 0; j < gridsize_tfp; j++)
			{
				prob_y[i * gridsize_tfp + j] = prob_z[i * gridsize_tfp + j];
			}
			for (int k = 0; k < gridsize_b; k++)
			{
				grid_y_d[i * gridsize_b + k] = CalcAutarkyOutput(grid_y_nd[i], grid_b[k], d_y, d_yy, d_b, d_bb, d_yb);
				grid_uD_yb[i * gridsize_b + k] = POWERFUN(max(clb, grid_y_d[i * gridsize_b + k]), 1 - rra) / (1 - rra);
			}
		}

		delete[] z;
		delete[] prob_z;

	}
	else
	{

		delete[] z;
		delete[] prob_z;

		return(EXIT_FAILURE);
	}

	return(EXIT_SUCCESS);

}

void parms_bsl_mod::matlab_display_cpu_parms()
{
	// display in matlab
	mexPrintf("Parameters:\n");
	mexPrintf("rra: %g\n", rra);
	mexPrintf("beta: %g\n", beta);
	mexPrintf("ggq_mDval: %g\n", ggq_mDval);
	mexPrintf("ggq_mlb: %g\n", ggq_mlb);
	mexPrintf("ggq_mub: %g\n", ggq_mub);
	mexPrintf("ggq_mmean: %g\n", ggq_mmean);
	mexPrintf("ggq_mstd: %g\n", ggq_mstd);
	mexPrintf("ggq_mmass: %g\n", ggq_mmass);
	mexPrintf("errtol_bisect: %g\n", errtol_bisect);
	mexPrintf("maxiters_bisect: %d\n", maxiters_bisect);
	mexPrintf("rH: %g\n", rH);
	mexPrintf("rL: %g\n", rL);
	mexPrintf("holding_cost: %g\n", holding_cost);
	mexPrintf("prob_meet_mkt_maker_ND: %g\n", prob_meet_mkt_maker_ND);
	mexPrintf("prob_meet_mkt_maker_D: %g\n", prob_meet_mkt_maker_D);
	mexPrintf("bargain_power_mkt_maker_ND: %g\n", bargain_power_mkt_maker_ND);
	mexPrintf("bargain_power_mkt_maker_D: %g\n", bargain_power_mkt_maker_D);
	mexPrintf("reenterprob: %g\n", reenterprob);
	mexPrintf("debt_z: %g\n", debt_z);
	mexPrintf("debt_m: %g\n", debt_m);
	mexPrintf("debt_writedown_reentry: %g\n", debt_writedown_reentry);
	mexPrintf("max_def_prob: %g\n", max_def_prob);
	mexPrintf("gridsize_tfp: %d\n", gridsize_tfp);
	mexPrintf("gridsize_b: %d\n", gridsize_b);
	mexPrintf("bmax: %g\n", bmax);
	mexPrintf("tfp_lb: %g\n", tfp_lb);
	mexPrintf("tfp_ub: %g\n", tfp_ub);
	mexPrintf("conv_chk_tfp_lb: %g\n", conv_chk_tfp_lb);
	mexPrintf("conv_chk_tfp_ub: %g\n", conv_chk_tfp_ub);
	mexPrintf("conv_chk_b_ub: %g\n", conv_chk_b_ub);
	mexPrintf("rho_tfp: %g\n", rho_tfp);
	mexPrintf("sigma_tfp: %g\n", sigma_tfp);
	mexPrintf("d_y: %g\n", d_y);
	mexPrintf("d_yy: %g\n", d_yy);
	mexPrintf("d_b: %g\n", d_b);
	mexPrintf("d_bb: %g\n", d_bb);
	mexPrintf("d_yb: %g\n", d_yb);
	mexPrintf("maxiters: %d\n", maxiters);
	mexPrintf("errtolV: %g\n", errtolV);
	mexPrintf("errtolQ: %g\n", errtolQ);
	mexPrintf("NumUpdateRegimes: %d\n", NumUpdateRegimes);
	for (int i = 0; i < NumUpdateRegimes - 1; i++)
	{
		mexPrintf("UpdateRegimeIter[%d]: %d\n", i, UpdateRegimeIter[i]);
	}
	for (int i = 0; i < NumUpdateRegimes; i++ )
	{
		mexPrintf("UpdateWeightOldV[%d]: %g\n", i, UpdateWeightOldV[i]);
		mexPrintf("UpdateWeightOldQ[%d]: %g\n", i, UpdateWeightOldQ[i]);
		mexPrintf("UpdateWeightOldDefProb[%d]: %g\n", i, UpdateWeightOldDefProb[i]);
	}
}