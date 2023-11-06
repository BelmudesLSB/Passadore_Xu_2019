#ifndef NORMALDIST_MEX_H
#define NORMALDIST_MEX_H

#include "realtype.h"

// REAL_TYPE normpdf_host(REAL_TYPE x, REAL_TYPE sigma);
// REAL_TYPE normcdf_host(REAL_TYPE x, REAL_TYPE sigma);

__host__ __device__ REAL_TYPE normpdf(REAL_TYPE x, REAL_TYPE sigma);
__host__ __device__ REAL_TYPE normcdf(REAL_TYPE x, REAL_TYPE sigma);

#endif // !NORMALDIST_H