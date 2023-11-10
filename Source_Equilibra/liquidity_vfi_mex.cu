#include "realtype.h"
#include "liquidity_vfi_mex.h"
#include "liquidity_mex_defs.h" // includes device constant definitions
#include "normaldist_mex.h"
#include <cmath>
#include<thrust/device_vector.h>
#include<thrust/device_ptr.h>
#include<thrust/extrema.h>
#include <cstdlib>
#include <math_constants.h>
#include "mex.h"

// GPU device constants: preferences
__constant__ REAL_TYPE RRA;
__constant__ REAL_TYPE BETA;
__constant__ REAL_TYPE ONE_MINUS_RRA; // 1-RRA
__constant__ REAL_TYPE U_SCALING; // 1/(1-RRA)

// GPU device constants: numerical parameters
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
__constant__ int CONV_CHK_B_UB_IDX;

// GPU device constants: liquidity discounts
__constant__ REAL_TYPE RH;
__constant__ REAL_TYPE RL;
__constant__ REAL_TYPE HOLDING_COST;
__constant__ REAL_TYPE PROB_HL;
__constant__ REAL_TYPE PROB_LH_ND;
__constant__ REAL_TYPE PROB_LH_D;

// GPU device constants: debt market
__constant__ REAL_TYPE REENTERPROB;
__constant__ REAL_TYPE DEBT_Z;
__constant__ REAL_TYPE DEBT_M;

// GPU device constants: grid dimensions
__constant__ int GRIDSIZE_Y;
__constant__ int GRIDSIZE_B;

// GPU device constants: grids
__constant__ REAL_TYPE GRID_Y_ND[MAX_GRIDSIZE_Y];
__constant__ REAL_TYPE GRID_B[MAX_GRIDSIZE_B];

// GPU device constants: interpolation
__constant__ int GRID_B_REENTER_IDX[MAX_GRIDSIZE_B];
__constant__ REAL_TYPE GRID_RECOVERY_FRAC[MAX_GRIDSIZE_B];

__device__ void ggq_topdown(REAL_TYPE& CENDupdate, REAL_TYPE& qHupdate, REAL_TYPE& qLupdate, REAL_TYPE& def_prob_update, REAL_TYPE& mdeftresh,
	int idx_y, int idx_b, REAL_TYPE& VD, REAL_TYPE* d_qH_ND, REAL_TYPE* d_qL_ND, REAL_TYPE* d_qH_D, REAL_TYPE* d_qL_D, REAL_TYPE* d_CVND, REAL_TYPE* d_def_prob)
{		
	// initialize at mub
	REAL_TYPE Vnow = CUDART_NEG_INF;
	REAL_TYPE Cnow, Wnow;
	REAL_TYPE Ccandidate, Wcandidate, Vcandidate;
	int i, ichoice;
	REAL_TYPE mnow = GGQ_MUB;

	// initialize at m1
	for (i = 0; i < GRIDSIZE_B; i++)
	{
		// cs0 = y - b*[m + (1 - m)*z] + qCH(y, b')*[b' - (1 - m)*b]
		// NOTE: we require cs0+mlb to be strictly positive!

		// ! Define consumption given state variables for a given x.

		Ccandidate = GRID_Y_ND[idx_y] - GRID_B[idx_b] * (DEBT_M + (1 - DEBT_M) * DEBT_Z) + d_qH_ND[idx_y * GRIDSIZE_B + i] * (GRID_B[i] - (1 - DEBT_M) * GRID_B[idx_b]);
		
		if (Ccandidate + mnow >= C_LB && (d_def_prob[idx_y * GRIDSIZE_B + i] <= MAX_DEF_PROB || GRID_B[i] <= (1 - DEBT_M) * GRID_B[idx_b]))
		{
			Wcandidate = BETA * d_CVND[idx_y * GRIDSIZE_B + i];
			Vcandidate = U_SCALING * POWERFUN(Ccandidate + mnow, ONE_MINUS_RRA) + Wcandidate; //mnow=GGQ_MLB
			// U_SCALING = 1/(1-RRA)
			if (Vcandidate > Vnow)
			{
				Cnow = Ccandidate;
				Wnow = Wcandidate;
				Vnow = Vcandidate;
				// record choices
				ichoice = i;
			}
			else if (Vcandidate == Vnow)
			{
				if (Ccandidate < Cnow)//choose lowest c
				{
					Cnow = Ccandidate;
					Wnow = Wcandidate;
					Vnow = Vcandidate;
					// record choices
					ichoice = i;
				}
			}
		}
	} // initialize
	if (!isfinite(Vnow))
	{
		// no feasible choice; always default
		CENDupdate = VD;
		def_prob_update = 1.0;	
		qHupdate = (1 - PROB_HL) * d_qH_D[idx_y * GRIDSIZE_B + idx_b] + PROB_HL * d_qL_D[idx_y * GRIDSIZE_B + idx_b];
		qLupdate = max(0.0, (1 - PROB_LH_D) * d_qL_D[idx_y * GRIDSIZE_B + idx_b] + PROB_LH_D * d_qH_D[idx_y * GRIDSIZE_B + idx_b] - HOLDING_COST);
		mdeftresh = GGQ_MUB;
	}
	else // ! At least we can try to search for candidates to replace (x_1,m_1)
	{

		int ichoice_prev = ichoice;
		REAL_TYPE M1;
		REAL_TYPE Vnow_M1;
		REAL_TYPE mnext, Cnext, Wnext;
		bool bContinue = true;

		CENDupdate = 0.0, qHupdate = 0.0, qLupdate = 0.0;

		if (VD >= Vnow)
		{
			// always default
			mdeftresh = GGQ_MUB;
		} else {
			mdeftresh = CUDART_NEG_INF; // infinite number will be used as a flag
		}

		while (bContinue)
		{

			M1 = max(GGQ_MLB, C_LB - Cnow); //! Lowest possible value of m to consider given current x_now.
			Vnow_M1 = U_SCALING * POWERFUN(Cnow + M1, ONE_MINUS_RRA) + Wnow; // ! Lowest possible value given the policy.
			mnext = M1;
			Cnext = Cnow;
			bContinue = false; // ! If there is no policy we still need to check whether default is preffered.

			// go through choices
			for (i = 0; i < GRIDSIZE_B; i++)
			{	// ! it does not matter if the choice of debt is only [m_2,m_1] we need to guarantee that d_def_prob[idx_y * GRIDSIZE_B + i] <= MAX_DEF_PROB.
				if (i != ichoice_prev && (d_def_prob[idx_y * GRIDSIZE_B + i] <= MAX_DEF_PROB || GRID_B[i] <= (1 - DEBT_M) * GRID_B[idx_b]))
				{
					Ccandidate = GRID_Y_ND[idx_y] - GRID_B[idx_b] * (DEBT_M + (1 - DEBT_M) * DEBT_Z) + d_qH_ND[idx_y * GRIDSIZE_B + i] * (GRID_B[i] - (1 - DEBT_M) * GRID_B[idx_b]);
					if (Ccandidate > Cnow) // only larger c is possibly a better choice, feasibility is automatically satisfied!
					{
						Wcandidate = BETA * d_CVND[idx_y * GRIDSIZE_B + i];
						Vcandidate = U_SCALING * POWERFUN(Ccandidate + M1, ONE_MINUS_RRA) + Wcandidate;
						if (Vcandidate > Vnow_M1) // ! check at the lowest. 
						{
							// ! Go over all the possible options, checking until we find the maximum m such that indifference.
							// find bisection point, save a register-- store in Vcandidate
							Vcandidate = bisect_zero(M1, mnow, Cnow, Ccandidate, (Wnow - Wcandidate) / U_SCALING);
							if (Vcandidate > mnext)
							{
								// NOTE: if the difference is too small, the root from the bisection may be equal to MUB. This will produce an error.
								ichoice = i;
								mnext = Vcandidate;
								Cnext = Ccandidate;
								Wnext = Wcandidate;
								bContinue = true;
							}
							else if (Vcandidate == mnext && Ccandidate < Cnext) // choose lowest consumption
							{
								ichoice = i;
								mnext = Vcandidate;
								Cnext = Ccandidate;
								Wnext = Wcandidate;
								bContinue = true;
							}
						} // is a better choice
					} // possible better choice
				}
			} // note, it's possible for this loop to have no feasible choice (e.g. q=0)
			
			if (VD >= U_SCALING * POWERFUN(Cnow + mnow, ONE_MINUS_RRA) + Wnow) // ! always default in region
			{
				
				CENDupdate += VD * (normcdf(mnow - GGQ_MMEAN, GGQ_MSTD) - normcdf(GGQ_MLB - GGQ_MMEAN, GGQ_MSTD));
				qHupdate += (normcdf(mnow - GGQ_MMEAN, GGQ_MSTD) - normcdf(GGQ_MLB - GGQ_MMEAN, GGQ_MSTD)) * ((1 - PROB_HL) * d_qH_D[idx_y * GRIDSIZE_B + idx_b] + PROB_HL * d_qL_D[idx_y * GRIDSIZE_B + idx_b]);
				qLupdate += (normcdf(mnow - GGQ_MMEAN, GGQ_MSTD) - normcdf(GGQ_MLB - GGQ_MMEAN, GGQ_MSTD)) * ((1 - PROB_LH_D) * d_qL_D[idx_y * GRIDSIZE_B + idx_b] + PROB_LH_D * d_qH_D[idx_y * GRIDSIZE_B + idx_b]);
				bContinue = false;
			}
			else if (VD <= U_SCALING * POWERFUN(Cnow + mnext, ONE_MINUS_RRA) + Wnow) // ! It might be the case that mnext is M1.
			{
				// ! never default in region until policy switches. (mnext is not M1)
				if (bContinue)
				{
					CENDupdate += gauss_legendre_CENDupdate(mnext, mnow, Cnow, Wnow);
					qHupdate += (normcdf(mnow - GGQ_MMEAN, GGQ_MSTD) - normcdf(mnext - GGQ_MMEAN, GGQ_MSTD)) *	(DEBT_M + (1 - DEBT_M) * DEBT_Z	+ (1 - DEBT_M) * ((1 - PROB_HL) * d_qH_ND[idx_y * GRIDSIZE_B + ichoice_prev] + PROB_HL * d_qL_ND[idx_y * GRIDSIZE_B + ichoice_prev]));
					qLupdate += (normcdf(mnow - GGQ_MMEAN, GGQ_MSTD) - normcdf(mnext - GGQ_MMEAN, GGQ_MSTD)) *	(DEBT_M + (1 - DEBT_M) * DEBT_Z	+ (1 - DEBT_M) * ((1 - PROB_LH_ND) * d_qL_ND[idx_y * GRIDSIZE_B + ichoice_prev] + PROB_LH_ND * d_qH_ND[idx_y * GRIDSIZE_B + ichoice_prev]));
					// update and continue algorithm
					ichoice_prev = ichoice;
					mnow = mnext;
					Cnow = Cnext;
					Wnow = Wnext;
				}
				else // ! (mnext is M1)
				{
					if (C_LB - Cnow > GGQ_MLB)
					{
						// (Cnow,Wnow) applies until M1=C_LB-Cnow>MLB
						CENDupdate += gauss_legendre_CENDupdate(M1, mnow, Cnow, Wnow);
						qHupdate += (normcdf(mnow - GGQ_MMEAN, GGQ_MSTD) - normcdf(M1 - GGQ_MMEAN, GGQ_MSTD)) *	(DEBT_M + (1 - DEBT_M) * DEBT_Z	+ (1 - DEBT_M) * ((1 - PROB_HL) * d_qH_ND[idx_y * GRIDSIZE_B + ichoice_prev] + PROB_HL * d_qL_ND[idx_y * GRIDSIZE_B + ichoice_prev]));
						qLupdate += (normcdf(mnow - GGQ_MMEAN, GGQ_MSTD) - normcdf(M1 - GGQ_MMEAN, GGQ_MSTD)) *	(DEBT_M + (1 - DEBT_M) * DEBT_Z	+ (1 - DEBT_M) * ((1 - PROB_LH_ND) * d_qL_ND[idx_y * GRIDSIZE_B + ichoice_prev] + PROB_LH_ND * d_qH_ND[idx_y * GRIDSIZE_B + ichoice_prev]));
						// reinitialize algorithm at the maximal value at M1
						bContinue = false;
						mnow = M1;
						Vnow = Vnow_M1;
						for (i = 0; i < GRIDSIZE_B; i++)
						{
							if (i != ichoice &&
								(d_def_prob[idx_y * GRIDSIZE_B + i] <= MAX_DEF_PROB || GRID_B[i] <= (1 - DEBT_M) * GRID_B[idx_b]))
							{
								Ccandidate = GRID_Y_ND[idx_y] - GRID_B[idx_b] * (DEBT_M + (1 - DEBT_M) * DEBT_Z)
									+ d_qH_ND[idx_y * GRIDSIZE_B + i] * (GRID_B[i] - (1 - DEBT_M) * GRID_B[idx_b]);
								if (Ccandidate + M1 > C_LB)
								{
									Wcandidate = BETA * d_CVND[idx_y * GRIDSIZE_B + i];
									Vcandidate = U_SCALING * POWERFUN(Ccandidate + M1, ONE_MINUS_RRA) + Wcandidate;
									if (Vcandidate > Vnow)
									{
										Cnow = Ccandidate;
										Wnow = Wcandidate;
										Vnow = Vcandidate;
										// record choices
										ichoice_prev = i;
										bContinue = true;
									}
									else if (Vcandidate == Vnow)
									{
										if (Ccandidate < Cnow)
										{
											Cnow = Ccandidate;
											Wnow = Wcandidate;
											Vnow = Vcandidate;
											// record choices
											ichoice_prev = i;
											bContinue = true;
										}
									}
								}
							}
						}
						if (!bContinue)
						{
							// all choices are not feasible, VD applies for m<=M1
							CENDupdate += VD * (normcdf(M1 - GGQ_MMEAN, GGQ_MSTD) - normcdf(GGQ_MLB - GGQ_MMEAN, GGQ_MSTD));
							qHupdate += (normcdf(M1 - GGQ_MMEAN, GGQ_MSTD) - normcdf(GGQ_MLB - GGQ_MMEAN, GGQ_MSTD)) *	((1 - PROB_HL) * d_qH_D[idx_y * GRIDSIZE_B + idx_b] + PROB_HL * d_qL_D[idx_y * GRIDSIZE_B + idx_b]);
							qLupdate += (normcdf(M1 - GGQ_MMEAN, GGQ_MSTD) - normcdf(GGQ_MLB - GGQ_MMEAN, GGQ_MSTD)) *	((1 - PROB_LH_D) * d_qL_D[idx_y * GRIDSIZE_B + idx_b] + PROB_LH_D * d_qH_D[idx_y * GRIDSIZE_B + idx_b]);
							mdeftresh = M1;
						}
					}
					else
					{
						// (Cnow,Wnow) applies until MLB
						CENDupdate += gauss_legendre_CENDupdate(GGQ_MLB, mnow, Cnow, Wnow);
						qHupdate += (normcdf(mnow - GGQ_MMEAN, GGQ_MSTD) - normcdf(GGQ_MLB - GGQ_MMEAN, GGQ_MSTD)) * (DEBT_M + (1 - DEBT_M) * DEBT_Z + (1 - DEBT_M) * ((1 - PROB_HL) * d_qH_ND[idx_y * GRIDSIZE_B + ichoice_prev] + PROB_HL * d_qL_ND[idx_y * GRIDSIZE_B + ichoice_prev]));
						qLupdate += (normcdf(mnow - GGQ_MMEAN, GGQ_MSTD) - normcdf(GGQ_MLB - GGQ_MMEAN, GGQ_MSTD)) * (DEBT_M + (1 - DEBT_M) * DEBT_Z + (1 - DEBT_M) * ((1 - PROB_LH_ND) * d_qL_ND[idx_y * GRIDSIZE_B + ichoice_prev] + PROB_LH_ND * d_qH_ND[idx_y * GRIDSIZE_B + ichoice_prev]));
						mdeftresh = GGQ_MLB;
					}
				}
			}
			else
			{
				// intermediate default
				Vcandidate = POWERFUN(ONE_MINUS_RRA * (VD - Wnow), U_SCALING) - Cnow; // def threshold
				CENDupdate += gauss_legendre_CENDupdate(Vcandidate, mnow, Cnow, Wnow) + VD * (normcdf(Vcandidate - GGQ_MMEAN, GGQ_MSTD) - normcdf(GGQ_MLB - GGQ_MMEAN, GGQ_MSTD));
				qHupdate += (normcdf(Vcandidate - GGQ_MMEAN, GGQ_MSTD) - normcdf(GGQ_MLB - GGQ_MMEAN, GGQ_MSTD)) *	((1 - PROB_HL) * d_qH_D[idx_y * GRIDSIZE_B + idx_b] + PROB_HL * d_qL_D[idx_y * GRIDSIZE_B + idx_b])	+ (normcdf(mnow - GGQ_MMEAN, GGQ_MSTD) - normcdf(Vcandidate - GGQ_MMEAN, GGQ_MSTD)) * (DEBT_M + (1 - DEBT_M) * DEBT_Z + (1 - DEBT_M) * ((1 - PROB_HL) * d_qH_ND[idx_y * GRIDSIZE_B + ichoice_prev] + PROB_HL * d_qL_ND[idx_y * GRIDSIZE_B + ichoice_prev]));
				qLupdate += (normcdf(Vcandidate - GGQ_MMEAN, GGQ_MSTD) - normcdf(GGQ_MLB - GGQ_MMEAN, GGQ_MSTD)) * 	((1 - PROB_LH_D) * d_qL_D[idx_y * GRIDSIZE_B + idx_b] + PROB_LH_D * d_qH_D[idx_y * GRIDSIZE_B + idx_b])	+ (normcdf(mnow - GGQ_MMEAN, GGQ_MSTD) - normcdf(Vcandidate - GGQ_MMEAN, GGQ_MSTD)) * (DEBT_M + (1 - DEBT_M) * DEBT_Z + (1 - DEBT_M) * ((1 - PROB_LH_ND) * d_qL_ND[idx_y * GRIDSIZE_B + ichoice_prev] + PROB_LH_ND * d_qH_ND[idx_y * GRIDSIZE_B + ichoice_prev]));
				bContinue = false;
				mdeftresh = Vcandidate;
			}

		} // while loop

		CENDupdate = CENDupdate / GGQ_MMASS;
		qHupdate = max(0.0, qHupdate) / GGQ_MMASS;
		def_prob_update = (normcdf(mdeftresh - GGQ_MMEAN, GGQ_MSTD) - normcdf(GGQ_MLB - GGQ_MMEAN, GGQ_MSTD)) / GGQ_MMASS;
		qLupdate = max(0.0, qLupdate - HOLDING_COST) / GGQ_MMASS;

	}

}

__device__ REAL_TYPE bisect_zero(REAL_TYPE a, REAL_TYPE b, REAL_TYPE c1, REAL_TYPE c2, REAL_TYPE W1_minus_W2)
{
	REAL_TYPE xlb = a, xub = b, xmid, fmid;
	REAL_TYPE flb = m_root_fun(xlb, c1, c2, W1_minus_W2);
	int i;
	for (i = 0; i < MAXITERS_BISECT; i++)
	{
		xmid = 0.5 * (xlb + xub);
		fmid = m_root_fun(xmid, c1, c2, W1_minus_W2);
		if (signbit(flb) == signbit(fmid))
		{
			// same sign
			xlb = xmid;
		}
		else
		{
			// different signs
			xub = xmid;
		}
		if (xub - xlb <= ERRTOL_BISECT)
		{
			break;
		}
	}

	return 0.5 * (xlb + xub);
}

__device__ REAL_TYPE m_root_fun(REAL_TYPE m, REAL_TYPE U1, REAL_TYPE U2, REAL_TYPE W1_minus_W2)
{
	return POWERFUN(U1 + m, ONE_MINUS_RRA) - POWERFUN(U2 + m, ONE_MINUS_RRA) + W1_minus_W2;
}

__device__ REAL_TYPE CENDupdatefun(REAL_TYPE m, REAL_TYPE U, REAL_TYPE W)
{
	return normpdf(m - GGQ_MMEAN, GGQ_MSTD)	* (U_SCALING * POWERFUN(U + m, ONE_MINUS_RRA) + W);
}

__device__ REAL_TYPE gauss_legendre_CENDupdate(REAL_TYPE m1, REAL_TYPE m2, REAL_TYPE U, REAL_TYPE W)
{
	REAL_TYPE alpha = 0.5 * (m1 + m2);
	REAL_TYPE beta = 0.5 * (m2 - m1);
	
	// 21 points
	return 	beta * (0.1460811336496904 * CENDupdatefun(alpha + beta * 0.0000000000000000, U, W)
		+ 0.1445244039899700 * CENDupdatefun(alpha + beta * -0.1455618541608951, U, W)
		+ 0.1445244039899700 * CENDupdatefun(alpha + beta * 0.1455618541608951, U, W)
		+ 0.1398873947910731 * CENDupdatefun(alpha + beta * -0.2880213168024011, U, W)
		+ 0.1398873947910731 * CENDupdatefun(alpha + beta * 0.2880213168024011, U, W)
		+ 0.1322689386333375 * CENDupdatefun(alpha + beta * -0.4243421202074388, U, W)
		+ 0.1322689386333375 * CENDupdatefun(alpha + beta * 0.4243421202074388, U, W)
		+ 0.1218314160537285 * CENDupdatefun(alpha + beta * -0.5516188358872198, U, W)
		+ 0.1218314160537285 * CENDupdatefun(alpha + beta * 0.5516188358872198, U, W)
		+ 0.1087972991671484 * CENDupdatefun(alpha + beta * -0.6671388041974123, U, W)
		+ 0.1087972991671484 * CENDupdatefun(alpha + beta * 0.6671388041974123, U, W)
		+ 0.0934444234560339 * CENDupdatefun(alpha + beta * -0.7684399634756779, U, W)
		+ 0.0934444234560339 * CENDupdatefun(alpha + beta * 0.7684399634756779, U, W)
		+ 0.0761001136283793 * CENDupdatefun(alpha + beta * -0.8533633645833173, U, W)
		+ 0.0761001136283793 * CENDupdatefun(alpha + beta * 0.8533633645833173, U, W)
		+ 0.0571344254268572 * CENDupdatefun(alpha + beta * -0.9200993341504008, U, W)
		+ 0.0571344254268572 * CENDupdatefun(alpha + beta * 0.9200993341504008, U, W)
		+ 0.0369537897708525 * CENDupdatefun(alpha + beta * -0.9672268385663063, U, W)
		+ 0.0369537897708525 * CENDupdatefun(alpha + beta * 0.9672268385663063, U, W)
		+ 0.0160172282577743 * CENDupdatefun(alpha + beta * -0.9937521706203895, U, W)
		+ 0.0160172282577743 * CENDupdatefun(alpha + beta * 0.9937521706203895, U, W));

}

/*
FORMATS:
1. VD(y,b), VND(y,b), q(y,b'), CVD(y,b), CVND(y,b')
*/

__global__ void vfi_iterate_policy(int *d_idx_bchoice, REAL_TYPE *d_qH_ND, REAL_TYPE *d_def_prob, REAL_TYPE *d_CVND)
// policy is for m = ggq_mmean, a value of -1 indicates an empty choice
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	volatile int idx_b = i % GRIDSIZE_B;
	volatile int idx_y = (i - idx_b) / GRIDSIZE_B;
	if (idx_y < GRIDSIZE_Y && idx_b < GRIDSIZE_B)
	{
		
		REAL_TYPE Vnow = CUDART_NEG_INF;
		REAL_TYPE Cnow;
		REAL_TYPE Ccandidate, Wcandidate, Vcandidate;
		int i, ichoice = -1;
		REAL_TYPE mnow = GGQ_MMEAN; // evaluate at mmean

		for (i = 0; i < GRIDSIZE_B; i++)
		{
			// cs0 = y - b*[m + (1 - m)*z] + qCH(y, b')*[b' - (1 - m)*b]
			// NOTE: we require cs0+mlb to be strictly positive!
			Ccandidate = GRID_Y_ND[idx_y] - GRID_B[idx_b] * (DEBT_M + (1 - DEBT_M) * DEBT_Z) + d_qH_ND[idx_y * GRIDSIZE_B + i] * (GRID_B[i] - (1 - DEBT_M) * GRID_B[idx_b]);
			if (Ccandidate + mnow >= C_LB && (d_def_prob[idx_y * GRIDSIZE_B + i] <= MAX_DEF_PROB || GRID_B[i] <= (1 - DEBT_M) * GRID_B[idx_b]))
			{
				Vcandidate = U_SCALING * POWERFUN(Ccandidate + mnow, ONE_MINUS_RRA) + BETA * d_CVND[idx_y * GRIDSIZE_B + i]; //mnow=GGQ_MLB
				// U_SCALING = 1/(1-RRA)
				if (Vcandidate > Vnow)
				{
					Cnow = Ccandidate;
					Vnow = Vcandidate;
					// record choices
					ichoice = i;
				}
				else if (Vcandidate == Vnow)
				{
					if (Ccandidate < Cnow)//choose lowest c
					{
						Cnow = Ccandidate;
						Vnow = Vcandidate;
						// record choices
						ichoice = i;
					}
				}
			}
		}

		d_idx_bchoice[idx_y * GRIDSIZE_B + idx_b] = ichoice;

	}
}

__global__ void vfi_iterate(REAL_TYPE *d_VD, REAL_TYPE *d_VNDupdate, REAL_TYPE *d_qH_ND_update, REAL_TYPE *d_qL_ND_update, 
	REAL_TYPE *d_defprob_update, REAL_TYPE *d_defthresh,
	REAL_TYPE *d_CVD, REAL_TYPE *d_qH_ND, REAL_TYPE *d_qL_ND, 
	REAL_TYPE *d_qH_D, REAL_TYPE *d_qL_D, REAL_TYPE *d_CVND, REAL_TYPE *d_uD_yb, REAL_TYPE *d_defprob)
// we don't store the debt policy here
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;

	volatile int idx_b = i % GRIDSIZE_B; // ! get index for bond.
	volatile int idx_y = (i - idx_b) / GRIDSIZE_B; // ! gets index for output.

	if (idx_y < GRIDSIZE_Y && idx_b < GRIDSIZE_B)
	{
		// update D problem
		REAL_TYPE VD = d_uD_yb[idx_y*GRIDSIZE_B + idx_b] + BETA*d_CVD[idx_y*GRIDSIZE_B + idx_b];
		REAL_TYPE CENDupdate, qHupdate, qLupdate, defprobupdate, defthresh;

		// partition and apply ggq algorithm on each sub-region.
		ggq_topdown(CENDupdate, qHupdate, qLupdate, defprobupdate, defthresh, idx_y, idx_b, VD,	d_qH_ND, d_qL_ND, d_qH_D, d_qL_D, d_CVND, d_defprob);

		d_VNDupdate[idx_y*GRIDSIZE_B + idx_b] = CENDupdate;
		d_qH_ND_update[idx_y*GRIDSIZE_B + idx_b] = qHupdate;
		d_qL_ND_update[idx_y*GRIDSIZE_B + idx_b] = qLupdate;
		d_defprob_update[idx_y*GRIDSIZE_B + idx_b] = defprobupdate;
		d_VD[idx_y*GRIDSIZE_B + idx_b] = VD;
		d_defthresh[idx_y*GRIDSIZE_B + idx_b] = defthresh;

	}
}

__global__ void vfi_interpolate(REAL_TYPE *ceCi, REAL_TYPE *qH_NDi, REAL_TYPE *qL_NDi, REAL_TYPE *ceC, REAL_TYPE *qH_ND, REAL_TYPE *qL_ND)
/*
interpolation based on writedown_rate:
ceD_haircut, qH_D, qL_D
interpolation based on writedown_rate and writedown_reenter:
ceC_haircut, qH_ND, qL_ND
*/
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	int idx_b = i % GRIDSIZE_B;
	int idx_y = (i - idx_b) / GRIDSIZE_B;

	if (i < GRIDSIZE_Y*GRIDSIZE_B)
	{
		ceCi[i] = ceC[idx_y*GRIDSIZE_B + GRID_B_REENTER_IDX[idx_b]];
		qH_NDi[i] = qH_ND[idx_y*GRIDSIZE_B + GRID_B_REENTER_IDX[idx_b]];
		qL_NDi[i] = qL_ND[idx_y*GRIDSIZE_B + GRID_B_REENTER_IDX[idx_b]];
	}
}

__global__ void vfi_update1(REAL_TYPE *d_prob_y, REAL_TYPE *d_EVD, REAL_TYPE *d_ceC, 
	REAL_TYPE *d_qH_ND, REAL_TYPE *d_qL_ND, REAL_TYPE *d_qH_D1, REAL_TYPE *d_qL_D1, REAL_TYPE *d_Edefprob,
	REAL_TYPE *d_VD, REAL_TYPE *d_VNDupdate, REAL_TYPE *d_qH_NDupdate, REAL_TYPE *d_qL_NDupdate, REAL_TYPE *d_defprob_update, 
	REAL_TYPE *d_qH_D0, REAL_TYPE *d_qL_D0)
	// note: no interpolation occurs here
	// ! we updated objects of the form (y',b'), we now take the expectations over y.
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	volatile int idx_b = i % GRIDSIZE_B;
	volatile int idx_y = (i - idx_b) / GRIDSIZE_B;
	if (idx_y < GRIDSIZE_Y && idx_b < GRIDSIZE_B)
	{
		REAL_TYPE tempf;

		// EVD
		tempf = 0;
		for (i = 0; i < GRIDSIZE_Y; i++)
		{
			tempf += d_prob_y[idx_y*GRIDSIZE_Y + i] * d_VD[i*GRIDSIZE_B + idx_b];
		}
		d_EVD[idx_y*GRIDSIZE_B + idx_b] = tempf;

		// ceC
		tempf = 0;
		for (i = 0; i < GRIDSIZE_Y; i++)
		{
			tempf += d_prob_y[idx_y*GRIDSIZE_Y + i] * d_VNDupdate[i*GRIDSIZE_B + idx_b];
		}
		d_ceC[idx_y*GRIDSIZE_B + idx_b] = tempf;

		// def prob
		tempf = 0;
		for (i = 0; i < GRIDSIZE_Y; i++)
		{
			tempf += d_prob_y[idx_y*GRIDSIZE_Y+i]* d_defprob_update[i*GRIDSIZE_B + idx_b];
		}
		d_Edefprob[idx_y*GRIDSIZE_B + idx_b] = tempf;

		// qND_H
		tempf = 0;
		for (i = 0; i < GRIDSIZE_Y; i++)
		{
			tempf += d_prob_y[idx_y*GRIDSIZE_Y+i]* d_qH_NDupdate[i*GRIDSIZE_B + idx_b];
		}
		d_qH_ND[idx_y*GRIDSIZE_B + idx_b] = tempf / (1 + RH);

		// qND_L
		tempf = 0;
		for (i = 0; i < GRIDSIZE_Y; i++)
		{
			tempf += d_prob_y[idx_y*GRIDSIZE_Y + i] * d_qL_NDupdate[i*GRIDSIZE_B + idx_b];
		}
		d_qL_ND[idx_y*GRIDSIZE_B + idx_b] = tempf / (1 + RL);

		// qD_H
		tempf = 0;
		for (i = 0; i < GRIDSIZE_Y; i++)
		{
			tempf += d_prob_y[idx_y*GRIDSIZE_Y + i] * ((1-PROB_HL)*d_qH_D0[i*GRIDSIZE_B + idx_b] + PROB_HL*d_qL_D0[i*GRIDSIZE_B + idx_b]);
		}
		d_qH_D1[idx_y*GRIDSIZE_B + idx_b] = tempf;

		// qD_L
		tempf = 0;
		for (i = 0; i < GRIDSIZE_Y; i++)
		{
			tempf += d_prob_y[idx_y*GRIDSIZE_Y + i] * max(0.0, (1 - PROB_LH_D)*d_qL_D0[i*GRIDSIZE_B + idx_b] + PROB_LH_D*d_qH_D0[i*GRIDSIZE_B + idx_b]  - HOLDING_COST);
		}
		d_qL_D1[idx_y*GRIDSIZE_B + idx_b] = tempf;

	}
}

__global__ void vfi_update2(REAL_TYPE *d_CVD, REAL_TYPE *d_CVDhaircut, REAL_TYPE *d_CVNDhaircut, REAL_TYPE *d_qH_D, REAL_TYPE *d_qL_D, REAL_TYPE *d_qH_Dhaircut, REAL_TYPE *d_qL_Dhaircut, REAL_TYPE *d_qH_NDhaircut, REAL_TYPE *d_qL_NDhaircut)
// update taking into account write downs
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < GRIDSIZE_Y*GRIDSIZE_B)
	{
		int idx_b = i % GRIDSIZE_B;
		d_CVD[i] = (1 - REENTERPROB)*d_CVDhaircut[i] + REENTERPROB*d_CVNDhaircut[i];
		d_qH_D[i] = max(0.0,(1 - REENTERPROB)*d_qH_Dhaircut[i] / (1 + RH) + GRID_RECOVERY_FRAC[idx_b] * REENTERPROB * d_qH_NDhaircut[i]);
		d_qL_D[i] = max(0.0,(1 - REENTERPROB)*d_qL_Dhaircut[i] / (1 + RL) + GRID_RECOVERY_FRAC[idx_b] * REENTERPROB * d_qL_NDhaircut[i]);
	}
}

__global__ void update_compute_errors(REAL_TYPE *d_x, REAL_TYPE *d_y, REAL_TYPE wx)
// new variable stored in x
// errors stored in y
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < GRIDSIZE_Y*GRIDSIZE_B)
	{
		REAL_TYPE x = d_x[i];
		REAL_TYPE y = d_y[i];
		d_x[i] = wx*x + (1 - wx)*y;
		int idx_b = i % GRIDSIZE_B;
		int idx_y = (i - idx_b)/GRIDSIZE_B;
		if (idx_y>=CONV_CHK_TFP_LB_IDX && idx_y<=CONV_CHK_TFP_UB_IDX && idx_b<= CONV_CHK_B_UB_IDX)
		{
			d_y[i] = abs(x - y);
		}
		else
		{
			d_y[i] = 0.0;
		}
	}
}

void vfi(parms_bsl_mod &p, REAL_TYPE wOldV, REAL_TYPE wOldQ, REAL_TYPE wOldDefProb, REAL_TYPE *d_prob_y, REAL_TYPE *d_uD_yb, 
	REAL_TYPE &err_CVD, REAL_TYPE &err_CVND, REAL_TYPE &err_qH_ND, REAL_TYPE &err_qL_ND, REAL_TYPE &err_qH_D, REAL_TYPE &err_qL_D, REAL_TYPE &err_defprob,
	REAL_TYPE *d_CVD, REAL_TYPE *d_CVND, REAL_TYPE *d_qH_ND, REAL_TYPE *d_qL_ND, REAL_TYPE *d_qH_D, REAL_TYPE *d_qL_D, REAL_TYPE *d_defprob, REAL_TYPE *d_VD, 
	REAL_TYPE *d_VNDupdate, REAL_TYPE *d_qH_ND_update, REAL_TYPE *d_qL_ND_update, REAL_TYPE *d_defprob_update, REAL_TYPE *d_defthresh,
	REAL_TYPE *d_EVD, REAL_TYPE *d_CVNDnew,
	REAL_TYPE *d_qH_NDnew, REAL_TYPE *d_qL_NDnew, REAL_TYPE *d_EqH_D, REAL_TYPE *d_EqL_D, REAL_TYPE *d_defprob_new,
	REAL_TYPE *d_CVNDi, REAL_TYPE *d_qH_NDi, REAL_TYPE *d_qL_NDi,
	REAL_TYPE *d_CVDnew, REAL_TYPE *d_qH_Dnew, REAL_TYPE *d_qL_Dnew)
	/*
	All vectors should have dimension GRIDSIZE_Y*GRIDSIZE_B
	// the device has plenty of memory
	*/
{	// ! Same number of blocks and threads per block in all kernels.

	int NUMTHREADS1 = 128;
	int NUMBLOCKS1 = (int) ceil((double)(p.gridsize_tfp*p.gridsize_b) / (double)NUMTHREADS1);

	// ! This runs the GGQ Algorithm:
	vfi_iterate <<< NUMBLOCKS1, NUMTHREADS1 >>> (d_VD, d_VNDupdate, d_qH_ND_update, d_qL_ND_update, d_defprob_update, d_defthresh, d_CVD, d_qH_ND, d_qL_ND,d_qH_D, d_qL_D, d_CVND, d_uD_yb, d_defprob); // inputs
	cudaDeviceSynchronize();

	// first round of updating
	int NUMTHREADS2 = 128;
	int NUMBLOCKS2 = (int)ceil((double)(p.gridsize_tfp*p.gridsize_b) / (double)NUMTHREADS2);

	// ! UPDATE Prices??
	vfi_update1 <<< NUMBLOCKS2, NUMTHREADS2 >>> (d_prob_y, d_EVD, d_CVNDnew, d_qH_NDnew, d_qL_NDnew, d_EqH_D, d_EqL_D, d_defprob_new, d_VD, d_VNDupdate, d_qH_ND_update, d_qL_ND_update, d_defprob_update, d_qH_D, d_qL_D);
	cudaDeviceSynchronize();

	// interpolate

	// ! Interpolation??
	vfi_interpolate <<< NUMBLOCKS2, NUMTHREADS2 >>> (d_CVNDi, d_qH_NDi, d_qL_NDi, d_CVNDnew, d_qH_NDnew, d_qL_NDnew);
	cudaDeviceSynchronize();

	// second round of updating

	// ! Renter probabilities ??
	vfi_update2 <<< NUMBLOCKS2, NUMTHREADS2 >>> (d_CVDnew, d_EVD, d_CVNDi, d_qH_Dnew, d_qL_Dnew, d_EqH_D, d_EqL_D, d_qH_NDi, d_qL_NDi);
	cudaDeviceSynchronize();

	// update and compute errors
	int NUMTHREADS3 = 128;
	int NUMBLOCKS3 = (int)ceil((double)(p.gridsize_tfp*p.gridsize_b) / (double)NUMTHREADS3);

	// ! Compute errors:
	update_compute_errors <<< NUMBLOCKS3, NUMTHREADS3 >>> (d_CVD, d_CVDnew, wOldV); // CVD
	update_compute_errors <<< NUMBLOCKS3, NUMTHREADS3 >>> (d_CVND, d_CVNDnew, wOldV); // CVND
	update_compute_errors <<< NUMBLOCKS3, NUMTHREADS3 >>> (d_defprob, d_defprob_new, wOldDefProb); // def prob
	update_compute_errors <<< NUMBLOCKS3, NUMTHREADS3 >>> (d_qH_ND, d_qH_NDnew, wOldQ); // qH_ND
	update_compute_errors <<< NUMBLOCKS3, NUMTHREADS3 >>> (d_qL_ND, d_qL_NDnew, wOldQ); // qL_ND
	update_compute_errors <<< NUMBLOCKS3, NUMTHREADS3 >>> (d_qH_D, d_qH_Dnew, wOldQ); // qH_D
	update_compute_errors <<< NUMBLOCKS3, NUMTHREADS3 >>> (d_qL_D, d_qL_Dnew, wOldQ); // qL_D

	cudaDeviceSynchronize();

	// reduce for norms
	thrust::device_ptr<REAL_TYPE> tempptr1 = thrust::device_pointer_cast(d_CVDnew);
	thrust::device_vector<REAL_TYPE>::iterator tempiter1 = thrust::max_element(tempptr1, tempptr1 + p.gridsize_tfp*p.gridsize_b);
	err_CVD = *tempiter1;

	thrust::device_ptr<REAL_TYPE> tempptr2 = thrust::device_pointer_cast(d_CVNDnew);
	thrust::device_vector<REAL_TYPE>::iterator tempiter2 = thrust::max_element(tempptr2, tempptr2 + p.gridsize_tfp*p.gridsize_b);
	err_CVND = *tempiter2;

	thrust::device_ptr<REAL_TYPE> tempptr3 = thrust::device_pointer_cast(d_qH_NDnew);
	thrust::device_vector<REAL_TYPE>::iterator tempiter3 = thrust::max_element(tempptr3, tempptr3 + p.gridsize_tfp*p.gridsize_b);
	err_qH_ND = *tempiter3;

	thrust::device_ptr<REAL_TYPE> tempptr4 = thrust::device_pointer_cast(d_qL_NDnew);
	thrust::device_vector<REAL_TYPE>::iterator tempiter4 = thrust::max_element(tempptr4, tempptr4 + p.gridsize_tfp*p.gridsize_b);
	err_qL_ND = *tempiter4;

	thrust::device_ptr<REAL_TYPE> tempptr5 = thrust::device_pointer_cast(d_qH_Dnew);
	thrust::device_vector<REAL_TYPE>::iterator tempiter5 = thrust::max_element(tempptr5, tempptr5 + p.gridsize_tfp*p.gridsize_b);
	err_qH_D = *tempiter5;

	thrust::device_ptr<REAL_TYPE> tempptr6 = thrust::device_pointer_cast(d_qL_Dnew);
	thrust::device_vector<REAL_TYPE>::iterator tempiter6 = thrust::max_element(tempptr6, tempptr6 + p.gridsize_tfp*p.gridsize_b);
	err_qL_D = *tempiter6;

	thrust::device_ptr<REAL_TYPE> tempptr7 = thrust::device_pointer_cast(d_defprob_new);
	thrust::device_vector<REAL_TYPE>::iterator tempiter7 = thrust::max_element(tempptr7, tempptr7 + p.gridsize_tfp*p.gridsize_b);
	err_defprob = *tempiter7;

	cudaDeviceSynchronize();
}

__global__ void initValueFuns(REAL_TYPE *d_CVD, REAL_TYPE *d_CVND, REAL_TYPE *d_qH_ND, REAL_TYPE *d_qL_ND, REAL_TYPE *d_qH_D, REAL_TYPE *d_qL_D,
	REAL_TYPE *d_defprob)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < GRIDSIZE_Y*GRIDSIZE_B)
	{
		d_CVD[i] = 0;
		d_CVND[i] = 0;
		d_qH_ND[i] = 0;
		d_qL_ND[i] = 0;
		d_qH_D[i] = 0;
		d_qL_D[i] = 0;
		d_defprob[i] = 0;
	}
}

void parms_bsl_mod::initDeviceConstants()
{
	REAL_TYPE tempf;
	// preferences
	cudaMemcpyToSymbol(RRA, &rra, 1 * sizeof(REAL_TYPE), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(BETA, &beta, 1 * sizeof(REAL_TYPE), 0, cudaMemcpyHostToDevice);
	tempf = 1 - rra;
	cudaMemcpyToSymbol(ONE_MINUS_RRA, &tempf, 1 * sizeof(REAL_TYPE), 0, cudaMemcpyHostToDevice);
	tempf = 1 / (1 - rra);
	cudaMemcpyToSymbol(U_SCALING, &tempf, 1 * sizeof(REAL_TYPE), 0, cudaMemcpyHostToDevice);
	// numerical parameters
	cudaMemcpyToSymbol(C_LB, &clb, 1 * sizeof(REAL_TYPE), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(GGQ_MLB, &ggq_mlb, 1 * sizeof(REAL_TYPE), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(GGQ_MUB, &ggq_mub, 1 * sizeof(REAL_TYPE), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(GGQ_MMEAN, &ggq_mmean, 1 * sizeof(REAL_TYPE), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(GGQ_MSTD, &ggq_mstd, 1 * sizeof(REAL_TYPE), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(GGQ_MMASS, &ggq_mmass, 1 * sizeof(REAL_TYPE), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(GGQ_MDEF, &ggq_mDval, 1 * sizeof(REAL_TYPE), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(ERRTOL_BISECT, &errtol_bisect, 1 * sizeof(REAL_TYPE), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(MAXITERS_BISECT, &maxiters_bisect, 1 * sizeof(int), 0, cudaMemcpyHostToDevice);
	// liquidity discounts
	cudaMemcpyToSymbol(RH, &rH, 1 * sizeof(REAL_TYPE), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(RL, &rL, 1 * sizeof(REAL_TYPE), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(HOLDING_COST, &holding_cost, 1 * sizeof(REAL_TYPE), 0, cudaMemcpyHostToDevice);
	// effective liquidity transitions
	cudaMemcpyToSymbol(PROB_HL, &prob_liqshock, 1 * sizeof(REAL_TYPE), 0, cudaMemcpyHostToDevice);
	tempf = prob_meet_mkt_maker_ND * (1 - bargain_power_mkt_maker_ND);
	cudaMemcpyToSymbol(PROB_LH_ND, &tempf, 1 * sizeof(REAL_TYPE), 0, cudaMemcpyHostToDevice);
	tempf = prob_meet_mkt_maker_D * (1 - bargain_power_mkt_maker_D);
	cudaMemcpyToSymbol(PROB_LH_D, &tempf, 1 * sizeof(REAL_TYPE), 0, cudaMemcpyHostToDevice);
	// debt market
	cudaMemcpyToSymbol(REENTERPROB, &reenterprob, 1 * sizeof(REAL_TYPE), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(DEBT_Z, &debt_z, 1 * sizeof(REAL_TYPE), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(DEBT_M, &debt_m, 1 * sizeof(REAL_TYPE), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(MAX_DEF_PROB, &max_def_prob, 1 * sizeof(REAL_TYPE), 0, cudaMemcpyHostToDevice);
	// grid dimensions
	cudaMemcpyToSymbol(GRIDSIZE_Y, &gridsize_tfp, 1 * sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(GRIDSIZE_B, &gridsize_b, 1 * sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(GRID_Y_ND, grid_y_nd, gridsize_tfp * sizeof(REAL_TYPE), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(GRID_B, grid_b, gridsize_b * sizeof(REAL_TYPE), 0, cudaMemcpyHostToDevice);

	// nearest point interpolation for reentry debt
	int* grid_b_reenter_idx = new int[gridsize_b];
	int tempidx;
	REAL_TYPE tempdiff;
	REAL_TYPE* grid_recovery_frac = new REAL_TYPE[gridsize_b];
	for (int idx_b = 0; idx_b < gridsize_b; idx_b++)
	{
		tempidx = 0;
		tempf = (idx_b > 0) ? grid_b[idx_b] * min(1 - debt_writedown_reentry, recovery_bmax / grid_b[idx_b]) : grid_b[idx_b] * (1 - debt_writedown_reentry);
		tempdiff = abs(grid_b[0] - tempf);
		for (int j = 1; j < gridsize_b; j++)
		{
			if (abs(grid_b[j] - tempf) < tempdiff)
			{
				tempdiff = abs(grid_b[j] - tempf);
				tempidx = j;
			}
		}
		grid_b_reenter_idx[idx_b] = tempidx;
		grid_recovery_frac[idx_b] = (idx_b > 0) ? min(1 - debt_writedown_reentry, recovery_bmax / grid_b[idx_b]) : (1 - debt_writedown_reentry);
	}

	cudaMemcpyToSymbol(GRID_B_REENTER_IDX, grid_b_reenter_idx, gridsize_b * sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(GRID_RECOVERY_FRAC, grid_recovery_frac, gridsize_b * sizeof(REAL_TYPE), 0, cudaMemcpyHostToDevice);
	delete[] grid_b_reenter_idx;
	delete[] grid_recovery_frac;

	// compute CONV_CHK_TFP_LB_IDX and CONV_CHK_TFP_UB_IDX;
	int conv_chk_tfp_lb_idx, conv_chk_tfp_ub_idx, conv_chk_b_ub_idx;
	for (int i = 0; i < gridsize_tfp; i++)
	{
		if (log(grid_y_nd[i]) >= conv_chk_tfp_lb)
		{
			conv_chk_tfp_lb_idx = i;
			break;
		}
	}
	for (int i = 0; i < gridsize_tfp; i++)
	{
		if (log(grid_y_nd[i]) <= conv_chk_tfp_ub)
		{
			conv_chk_tfp_ub_idx = i;
		}
		else
		{
			break;
		}
	}
	for (int i = 0; i < gridsize_b; i++)
	{
		if (grid_b[i] <= conv_chk_b_ub)
		{
			conv_chk_b_ub_idx = i;
		}
		else
		{
			break;
		}
	}

	cudaMemcpyToSymbol(CONV_CHK_TFP_LB_IDX, &conv_chk_tfp_lb_idx, 1 * sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(CONV_CHK_TFP_UB_IDX, &conv_chk_tfp_ub_idx, 1 * sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(CONV_CHK_B_UB_IDX, &conv_chk_b_ub_idx, 1 * sizeof(int), 0, cudaMemcpyHostToDevice);

	cudaDeviceSynchronize();

}

int SolveModel(REAL_TYPE* h_CVD, REAL_TYPE* h_CVND, REAL_TYPE* h_qH_ND, REAL_TYPE* h_qL_ND, REAL_TYPE* h_qH_D, REAL_TYPE*h_qL_D, 
	REAL_TYPE* h_defprob, REAL_TYPE* h_defthresh, int *h_idx_bchoice,  
	int &iter, REAL_TYPE &err_qH_ND, REAL_TYPE &err_qL_ND, REAL_TYPE &err_qH_D, REAL_TYPE &err_qL_D, REAL_TYPE &err_defprob, 
	REAL_TYPE &err_CVD, REAL_TYPE &err_CVND,
	parms_bsl_mod &p, int useDevice, bool bDisplayProgress, int DisplayInterval)
{

	if (bDisplayProgress) {
		p.matlab_display_cpu_parms();
	}

	if (p.checkparms() == EXIT_FAILURE)
	{
		mexPrintf("Error with parameter(s).");
		return(EXIT_FAILURE);
	}

	// initialize gpu
	int numDevices = 0;
	cudaGetDeviceCount(&numDevices);
	if (useDevice >= numDevices)
	{
		mexPrintf("Error: program requested device %d; only %d device(s) are available.\n", useDevice, numDevices);
		return(EXIT_FAILURE);
	} else {
		mexPrintf("Running with GPU device %d of %d.\n", useDevice, numDevices);
	}

	cudaSetDevice(useDevice);

	p.initDeviceConstants(); // ! Load constants to the GPU.

	// allocate device memory:
	REAL_TYPE *d_CVD, *d_CVND, *d_qH_ND, *d_qL_ND, *d_qH_D, *d_qL_D, *d_defprob;
	REAL_TYPE *d_VD, *d_VNDupdate, *d_qH_ND_update, *d_qL_ND_update, *d_EVD, *d_CVNDnew,
		*d_qH_NDnew, *d_qL_NDnew, *d_EqH_D, *d_EqL_D,
		*d_CVNDi, *d_qH_NDi, *d_qL_NDi,
		*d_CVDnew, *d_qH_Dnew, *d_qL_Dnew, *d_defprob_new, *d_defprob_update, *d_defthresh;
	REAL_TYPE *d_prob_y, *d_uD_yb;
	int* d_idx_bchoice;

	cudaMalloc((void**)&d_defprob, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE));
	cudaMalloc((void**)&d_defprob_new, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE));
	cudaMalloc((void**)&d_defprob_update, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE));
	cudaMalloc((void**)&d_defthresh, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE));
	cudaMalloc((void**)&d_CVD, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE));
	cudaMalloc((void**)&d_CVND, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE));
	cudaMalloc((void**)&d_qH_ND, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE));
	cudaMalloc((void**)&d_qL_ND, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE));
	cudaMalloc((void**)&d_qH_D, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE));
	cudaMalloc((void**)&d_qL_D, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE));
	cudaMalloc((void**)&d_VD, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE));
	cudaMalloc((void**)&d_VNDupdate, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE));
	cudaMalloc((void**)&d_qH_ND_update, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE));
	cudaMalloc((void**)&d_qL_ND_update, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE));
	cudaMalloc((void**)&d_EVD, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE));
	cudaMalloc((void**)&d_CVNDnew, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE));
	cudaMalloc((void**)&d_qH_NDnew, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE));
	cudaMalloc((void**)&d_qL_NDnew, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE));
	cudaMalloc((void**)&d_EqH_D, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE));
	cudaMalloc((void**)&d_EqL_D, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE));
	cudaMalloc((void**)&d_CVNDi, p.gridsize_tfp*p.gridsize_b * sizeof(REAL_TYPE));
	cudaMalloc((void**)&d_qH_NDi, p.gridsize_tfp*p.gridsize_b * sizeof(REAL_TYPE));
	cudaMalloc((void**)&d_qL_NDi, p.gridsize_tfp*p.gridsize_b * sizeof(REAL_TYPE));
	cudaMalloc((void**)&d_CVDnew, p.gridsize_tfp*p.gridsize_b * sizeof(REAL_TYPE));
	cudaMalloc((void**)&d_qH_Dnew, p.gridsize_tfp*p.gridsize_b * sizeof(REAL_TYPE));
	cudaMalloc((void**)&d_qL_Dnew, p.gridsize_tfp*p.gridsize_b * sizeof(REAL_TYPE));
	cudaMalloc((void**)&d_prob_y, p.gridsize_tfp*p.gridsize_tfp*sizeof(REAL_TYPE));
	cudaMemcpy(d_prob_y, p.prob_y, p.gridsize_tfp*p.gridsize_tfp*sizeof(REAL_TYPE), cudaMemcpyHostToDevice);
	cudaMalloc((void**)&d_uD_yb, p.gridsize_tfp*p.gridsize_b * sizeof(REAL_TYPE));
	cudaMemcpy(d_uD_yb, p.grid_uD_yb, p.gridsize_tfp*p.gridsize_b * sizeof(REAL_TYPE), cudaMemcpyHostToDevice);
	cudaMalloc((void**)&d_idx_bchoice, p.gridsize_tfp * p.gridsize_b * sizeof(int));

	cudaDeviceSynchronize();
	initValueFuns <<< (int)ceil((double)(p.gridsize_tfp*p.gridsize_b)/(double)128) , 128 >>> (d_CVD, d_CVND, d_qH_ND, d_qL_ND, d_qH_D, d_qL_D, d_defprob);
	cudaDeviceSynchronize();

	size_t free_cuda_mem, total_cuda_mem;
	cudaMemGetInfo(&free_cuda_mem, &total_cuda_mem);

	if (bDisplayProgress) {
		mexPrintf("Memory allocated for value function iteration. Total GPU memory: %g GBs. Free memory: %g GBs.\n",
			round(100 * (double)total_cuda_mem / 1073741824) / 100, round(100 * (double)free_cuda_mem / 1073741824) / 100);
	}


	cudaEvent_t start, stop;

	float milliseconds_taken;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start);

	iter = 0;
	int Updatek = 0;
	bool bConverged = false;
	REAL_TYPE err_q;
	REAL_TYPE wOldV, wOldQ, wOldDefProb;

	do
	{
		if (Updatek < p.NumUpdateRegimes - 1)
		{
			if (iter >= p.UpdateRegimeIter[Updatek])
			{
				Updatek += 1;
			}
		}

		iter += 1;
		wOldV = p.UpdateWeightOldV[Updatek];
		wOldQ = p.UpdateWeightOldQ[Updatek];
		wOldDefProb = p.UpdateWeightOldDefProb[Updatek];

		vfi(p, wOldV, wOldQ, wOldDefProb, d_prob_y, d_uD_yb, err_CVD, err_CVND, err_qH_ND, err_qL_ND, err_qH_D, err_qL_D, err_defprob,
			d_CVD, d_CVND, d_qH_ND, d_qL_ND, d_qH_D, d_qL_D, d_defprob,	d_VD, d_VNDupdate, d_qH_ND_update, d_qL_ND_update, d_defprob_update, d_defthresh,
			d_EVD, d_CVNDnew, d_qH_NDnew, d_qL_NDnew, d_EqH_D, d_EqL_D, d_defprob_new, d_CVNDi, d_qH_NDi, d_qL_NDi,	d_CVDnew, d_qH_Dnew, d_qL_Dnew);

		err_q = err_qH_ND;
		if (err_qL_ND > err_q)
		{
			err_q = err_qL_ND;
		}
		if (err_qH_D > err_q)
		{
			err_q = err_qH_D;
		}
		if (err_qL_D > err_q)
		{
			err_q = err_qL_D;
		}
		if (err_defprob > err_q)
		{
			err_q = err_defprob;
		}

		bConverged = err_q <= p.errtolQ && err_CVD <= p.errtolV && err_CVND <= p.errtolV;

		if (iter % DisplayInterval == 0 && bDisplayProgress)
		{
			cudaEventRecord(stop);
			cudaEventSynchronize(stop);
			cudaEventElapsedTime(&milliseconds_taken, start, stop);
			mexPrintf("iter %d, err(defprob)=%g, err(qs)=%g, err(CVD)=%g, err(CVND)=%g. %g seconds elapsed.\n", iter, err_defprob, err_q, err_CVD, err_CVND, milliseconds_taken / 1000.0);
		}

	} while (iter <= p.maxiters && !bConverged);

	// obtain policies at m = mmean // ! Why?

	cudaDeviceSynchronize();
	int NUMTHREADS1 = 128;
	int NUMBLOCKS1 = (int)ceil((double)(p.gridsize_tfp * p.gridsize_b) / (double)NUMTHREADS1);
	vfi_iterate_policy <<< NUMBLOCKS1, NUMTHREADS1 >>> (d_idx_bchoice, d_qH_ND, d_defprob, d_CVND);
	cudaDeviceSynchronize();

	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&milliseconds_taken, start, stop);

	if (bDisplayProgress) {
		mexPrintf("Program stopped after %d iterations. %g seconds elapsed.\n", iter, milliseconds_taken / 1000.0);
		mexPrintf("iter %d, err(defprob)=%g, err(qs)=%g, err(CVD)=%g, err(CVND)=%g. %g seconds elapsed.\n", iter, err_defprob, err_q, err_CVD, err_CVND, milliseconds_taken / 1000.0);
	}

	cudaMemcpy(h_defprob,d_defprob, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_CVD,d_CVD, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_CVND, d_CVND, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_qH_ND, d_qH_ND, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_qL_ND, d_qL_ND, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_qH_D, d_qH_D, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_qL_D, d_qL_D, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_defthresh, d_defthresh, p.gridsize_tfp*p.gridsize_b*sizeof(REAL_TYPE), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_idx_bchoice, d_idx_bchoice, p.gridsize_tfp * p.gridsize_b * sizeof(int), cudaMemcpyDeviceToHost);

	// free device memory
	cudaFree(d_defprob);
	cudaFree(d_defprob_update);
	cudaFree(d_defthresh);
	cudaFree(d_defprob_new);
	cudaFree(d_CVD);
	cudaFree(d_CVND);
	cudaFree(d_qH_ND);
	cudaFree(d_qL_ND);
	cudaFree(d_qH_D);
	cudaFree(d_qL_D);
	cudaFree(d_VD);
	cudaFree(d_VNDupdate);
	cudaFree(d_qH_ND_update);
	cudaFree(d_qL_ND_update);
	cudaFree(d_EVD);
	cudaFree(d_CVNDnew);
	cudaFree(d_qH_NDnew);
	cudaFree(d_qL_NDnew);
	cudaFree(d_EqH_D);
	cudaFree(d_EqL_D);
	cudaFree(d_CVNDi);
	cudaFree(d_qH_NDi);
	cudaFree(d_qL_NDi);
	cudaFree(d_CVDnew);
	cudaFree(d_qH_Dnew);
	cudaFree(d_qL_Dnew);
	cudaFree(d_prob_y);
	cudaFree(d_uD_yb);
	cudaFree(d_idx_bchoice);

	return(EXIT_SUCCESS);

}