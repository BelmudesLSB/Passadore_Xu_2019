#ifndef LIQUIDITY_OUTPUT_TO_MAT_H
#define LIQUIDITY_OUTPUT_TO_MAT_H

#include "mex.h"
#include "liquidity_mex_defs.h"

#ifndef NUM_PARM_FIELDS_TO_OUTPUT
// #define NUM_PARM_FIELDS_TO_OUTPUT 47 
#define NUM_PARM_FIELDS_TO_OUTPUT 2 
#endif


template <class T> void write_scalar_to_mat_struc(mxArray* pm, const char* fieldname, T write_val);
void save_parms_to_mat_struc(mxArray* pOut, parms_bsl_mod& p);
template <class T> void write_vector_to_mat_struc(mxArray* pm, const char* fieldname, T* write_val, int num_to_write);

#endif