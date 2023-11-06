#include "normaldist_mex.h"
#include "realtype.h"
#include<cmath>

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257390
#endif // M_2_SQRTPI

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif // !M_SQRT2

/*
REAL_TYPE normpdf_host(REAL_TYPE x, REAL_TYPE sigma) // mu=0, assumes that sigma>0
{
	// M_2_SQRTPI=2/sqrt(pi)=1.12837916709551257390
	return exp(-0.5 * x * x / (sigma * sigma)) * M_2_SQRTPI / (2 * sigma * M_SQRT2);
}

REAL_TYPE normcdf_host(REAL_TYPE x, REAL_TYPE sigma) // mu=0
{
	// c++'s errf(x) gives the integral of 2/sqrt(pi)*exp(-t^2) for 0<t<x
	return 0.5 * (1 + erf(x / (sigma * M_SQRT2)));
}
*/

__host__ __device__ REAL_TYPE normpdf(REAL_TYPE x, REAL_TYPE sigma) // mu=0, assumes that sigma>0
{
	// M_2_SQRTPI=2/sqrt(pi)=1.12837916709551257390
	return exp(-0.5*x*x/(sigma*sigma))*M_2_SQRTPI/(2*sigma*M_SQRT2);
}

__host__ __device__ REAL_TYPE normcdf(REAL_TYPE x, REAL_TYPE sigma) // mu=0
{
	// c++'s errf(x) gives the integral of 2/sqrt(pi)*exp(-t^2) for 0<t<x
	return 0.5*(1+erf(x/(sigma*M_SQRT2)));
}