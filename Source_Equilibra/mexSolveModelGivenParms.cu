#include "mex.h"
#include <cstdio>
#include <cstdlib>
#include "liquidity_mex_defs.h"
#include "liquidity_vfi_mex.h"
#include "matrix.h"


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
// see https://www.mathworks.com/matlabcentral/answers/94226-why-do-i-receive-the-error-mex-file-entry-point-is-missing (for mex files, not mexcuda)
// matlab syntax: output = mexSolveModelGivenParms(parm_struc, useDevice, DisplayInterval)


	if (nrhs != 3)
	{
		mexErrMsgIdAndTxt("mexSolveModelGivenParms:nrhs", "Correct format: output = mexSolveModelGivenParms(parm_struc, useDevice, DisplayInterval)\n");
	}

	if (nlhs != 1)
	{
		mexErrMsgIdAndTxt("mexSolveModelGivenParms:nlhs", "Correct format: output = mexSolveModelGivenParms(parm_struc, useDevice, DisplayInterval)\n");
	}

	parms_bsl_mod p;

	int useDevice = 0, DisplayInterval = 1;
	bool bDisplayProgress = true;

	if (!mxIsStruct(prhs[0])) {
		mexErrMsgIdAndTxt("mexSolveModelGivenParms:InputNotStruc", "Input[0] must be a structure.\n");
	}

	if (!mxIsDouble(prhs[1])) {
		mexErrMsgIdAndTxt("mexSolveModelGivenParms:InputNotDouble", "Input[1] must be a number.\n");
	}
	else {
		useDevice = (int)mxGetScalar(prhs[1]);
	}

	if (!mxIsDouble(prhs[2])) {
		mexErrMsgIdAndTxt("mexSolveModelGivenParms:InputNotDouble", "Input[2] must be a number.\n");
	}
	else {
		DisplayInterval = (int)mxGetScalar(prhs[2]);
		if (DisplayInterval >= 1) {
			bDisplayProgress = true;
		}
		else
		{
			bDisplayProgress = false;
		}
	}

	if (p.readparms(prhs[0]) == EXIT_SUCCESS)
	{

		// mexPrintf("Successfully read parameters.\n");
		// p.matlab_display_cpu_parms();

		// run vfi

		REAL_TYPE *h_CVD, *h_CVND, *h_qH_ND, *h_qL_ND, *h_qH_D, *h_qL_D, *h_defprob, *h_defthresh;
		int* h_idx_bchoice, iter;
		REAL_TYPE err_qH_ND, err_qL_ND, err_qH_D, err_qL_D, err_defprob, err_CVD, err_CVND;

		h_defprob = new REAL_TYPE[p.gridsize_tfp * p.gridsize_b];
		h_defthresh = new REAL_TYPE[p.gridsize_tfp * p.gridsize_b];
		h_CVD = new REAL_TYPE[p.gridsize_tfp * p.gridsize_b];
		h_CVND = new REAL_TYPE[p.gridsize_tfp * p.gridsize_b];
		h_qH_ND = new REAL_TYPE[p.gridsize_tfp * p.gridsize_b];
		h_qL_ND = new REAL_TYPE[p.gridsize_tfp * p.gridsize_b];
		h_qH_D = new REAL_TYPE[p.gridsize_tfp * p.gridsize_b];
		h_qL_D = new REAL_TYPE[p.gridsize_tfp * p.gridsize_b];
		h_idx_bchoice = new int[p.gridsize_tfp * p.gridsize_b];

		if (SolveModel(h_CVD, h_CVND, h_qH_ND, h_qL_ND, h_qH_D, h_qL_D, h_defprob, h_defthresh, h_idx_bchoice,
			iter, err_qH_ND, err_qL_ND, err_qH_D, err_qL_D, err_defprob, err_CVD, err_CVND,
			p, useDevice, bDisplayProgress, DisplayInterval) == EXIT_FAILURE) {
			mexErrMsgIdAndTxt("mexSolveModelGivenParms:SolveModel", "Error in SolveModel.\n");
		}

		// output to matlab
		const char* fnames[20] = {"iter", "err_qH_ND", "err_qL_ND", "err_qH_D", "err_qL_D", "err_defprob", "err_CVD", "err_CVND",
			"CVD","CVND","qH_ND","qL_ND","qH_D","qL_D","defprob","defthresh","idx_bchoice","prob_y","grid_y_nd","grid_y_d"};

		plhs[0] = mxCreateStructMatrix(1, 1, 20, fnames);

		mxArray *mxArray1;
		double *tmp_dbl;

		mxArray1 = mxCreateDoubleMatrix(1, 1, mxREAL);
		tmp_dbl = mxGetPr(mxArray1);
		tmp_dbl[0] = (double)iter;
		mxSetField(plhs[0], 0, "iter", mxArray1);

		mxArray1 = mxCreateDoubleMatrix(1, 1, mxREAL);
		tmp_dbl = mxGetPr(mxArray1);
		tmp_dbl[0] = (double)err_qH_ND;
		mxSetField(plhs[0], 0, "err_qH_ND", mxArray1);

		mxArray1 = mxCreateDoubleMatrix(1, 1, mxREAL);
		tmp_dbl = mxGetPr(mxArray1);
		tmp_dbl[0] = (double)err_qL_ND;
		mxSetField(plhs[0], 0, "err_qL_ND", mxArray1);

		mxArray1 = mxCreateDoubleMatrix(1, 1, mxREAL);
		tmp_dbl = mxGetPr(mxArray1);
		tmp_dbl[0] = (double)err_qH_D;
		mxSetField(plhs[0], 0, "err_qH_D", mxArray1);

		mxArray1 = mxCreateDoubleMatrix(1, 1, mxREAL);
		tmp_dbl = mxGetPr(mxArray1);
		tmp_dbl[0] = (double)err_qL_D;
		mxSetField(plhs[0], 0, "err_qL_D", mxArray1);

		mxArray1 = mxCreateDoubleMatrix(1, 1, mxREAL);
		tmp_dbl = mxGetPr(mxArray1);
		tmp_dbl[0] = (double)err_defprob;
		mxSetField(plhs[0], 0, "err_defprob", mxArray1);

		mxArray1 = mxCreateDoubleMatrix(1, 1, mxREAL);
		tmp_dbl = mxGetPr(mxArray1);
		tmp_dbl[0] = (double)err_CVD;
		mxSetField(plhs[0], 0, "err_CVD", mxArray1);

		mxArray1 = mxCreateDoubleMatrix(1, 1, mxREAL);
		tmp_dbl = mxGetPr(mxArray1);
		tmp_dbl[0] = (double)err_CVND;
		mxSetField(plhs[0], 0, "err_CVND", mxArray1);

		mxArray1 = mxCreateDoubleMatrix(p.gridsize_tfp * p.gridsize_b, 1, mxREAL);
		tmp_dbl = mxGetPr(mxArray1);
		for (int i = 0; i < p.gridsize_tfp * p.gridsize_b; i++)
		{
			tmp_dbl[i] = (double)h_CVD[i];
		}
		mxSetField(plhs[0], 0, "CVD", mxArray1);

		mxArray1 = mxCreateDoubleMatrix(p.gridsize_tfp * p.gridsize_b, 1, mxREAL);
		tmp_dbl = mxGetPr(mxArray1);
		for (int i = 0; i < p.gridsize_tfp * p.gridsize_b; i++)
		{
			tmp_dbl[i] = (double)h_CVND[i];
		}
		mxSetField(plhs[0], 0, "CVND", mxArray1);

		mxArray1 = mxCreateDoubleMatrix(p.gridsize_tfp * p.gridsize_b, 1, mxREAL);
		tmp_dbl = mxGetPr(mxArray1);
		for (int i = 0; i < p.gridsize_tfp * p.gridsize_b; i++)
		{
			tmp_dbl[i] = (double)h_qH_ND[i];
		}
		mxSetField(plhs[0], 0, "qH_ND", mxArray1);

		mxArray1 = mxCreateDoubleMatrix(p.gridsize_tfp * p.gridsize_b, 1, mxREAL);
		tmp_dbl = mxGetPr(mxArray1);
		for (int i = 0; i < p.gridsize_tfp * p.gridsize_b; i++)
		{
			tmp_dbl[i] = (double)h_qL_ND[i];
		}
		mxSetField(plhs[0], 0, "qL_ND", mxArray1);

		mxArray1 = mxCreateDoubleMatrix(p.gridsize_tfp * p.gridsize_b, 1, mxREAL);
		tmp_dbl = mxGetPr(mxArray1);
		for (int i = 0; i < p.gridsize_tfp * p.gridsize_b; i++)
		{
			tmp_dbl[i] = (double)h_qH_D[i];
		}
		mxSetField(plhs[0], 0, "qH_D", mxArray1);

		mxArray1 = mxCreateDoubleMatrix(p.gridsize_tfp * p.gridsize_b, 1, mxREAL);
		tmp_dbl = mxGetPr(mxArray1);
		for (int i = 0; i < p.gridsize_tfp * p.gridsize_b; i++)
		{
			tmp_dbl[i] = (double)h_qL_D[i];
		}
		mxSetField(plhs[0], 0, "qL_D", mxArray1);

		mxArray1 = mxCreateDoubleMatrix(p.gridsize_tfp * p.gridsize_b, 1, mxREAL);
		tmp_dbl = mxGetPr(mxArray1);
		for (int i = 0; i < p.gridsize_tfp * p.gridsize_b; i++)
		{
			tmp_dbl[i] = (double)h_defprob[i];
		}
		mxSetField(plhs[0], 0, "defprob", mxArray1);

		mxArray1 = mxCreateDoubleMatrix(p.gridsize_tfp * p.gridsize_b, 1, mxREAL);
		tmp_dbl = mxGetPr(mxArray1);
		for (int i = 0; i < p.gridsize_tfp * p.gridsize_b; i++)
		{
			tmp_dbl[i] = (double)h_defthresh[i];
		}
		mxSetField(plhs[0], 0, "defthresh", mxArray1);

		mxArray1 = mxCreateDoubleMatrix(p.gridsize_tfp * p.gridsize_b, 1, mxREAL);
		tmp_dbl = mxGetPr(mxArray1);
		for (int i = 0; i < p.gridsize_tfp * p.gridsize_b; i++)
		{
			tmp_dbl[i] = (double)h_idx_bchoice[i];
		}
		mxSetField(plhs[0], 0, "idx_bchoice", mxArray1);

		mxArray1 = mxCreateDoubleMatrix(p.gridsize_tfp * p.gridsize_tfp, 1, mxREAL);
		tmp_dbl = mxGetPr(mxArray1);
		for (int i = 0; i < p.gridsize_tfp * p.gridsize_tfp; i++)
		{
			tmp_dbl[i] = (double)p.prob_y[i];
		}
		mxSetField(plhs[0], 0, "prob_y", mxArray1);

		mxArray1 = mxCreateDoubleMatrix(p.gridsize_tfp, 1, mxREAL);
		tmp_dbl = mxGetPr(mxArray1);
		for (int i = 0; i < p.gridsize_tfp; i++)
		{
			tmp_dbl[i] = (double)p.grid_y_nd[i];
		}
		mxSetField(plhs[0], 0, "grid_y_nd", mxArray1);

		mxArray1 = mxCreateDoubleMatrix(p.gridsize_tfp*p.gridsize_b, 1, mxREAL);
		tmp_dbl = mxGetPr(mxArray1);
		for (int i = 0; i < p.gridsize_tfp*p.gridsize_b; i++)
		{
			tmp_dbl[i] = (double)p.grid_y_d[i];
		}
		mxSetField(plhs[0], 0, "grid_y_d", mxArray1);

	}
	else
	{
		mexErrMsgIdAndTxt("mexSolveModelGivenParms:ReadParms", "Error reading parameters.\n");
	}

}