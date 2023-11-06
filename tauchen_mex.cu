#include<cmath>
#include<cstdlib>
#include<cstdio>
#include "realtype.h"
#include "normaldist_mex.h"
#include "mex.h"

// #define DEBUG_IDENTITY_MATRIX

int tauchen(REAL_TYPE *x, REAL_TYPE *Px, const REAL_TYPE rho, const REAL_TYPE sigma, const int N, const REAL_TYPE xlb, const REAL_TYPE xub)
// Model:    x(t+1) = rho x(t) + sigma*e(t+1)
// evenly spaced grid between xlb and xub
// the space for x and Px should be already allocated
{

    if (abs(rho)<1 && sigma>0 && N>=2 && xlb<0.0 && xub>0.0)
    {

        REAL_TYPE dx = (xub - xlb)/((REAL_TYPE)N - 1.0);
        REAL_TYPE x1, x2, mean_x;
        x[0] = xlb;
        for (int i = 1; i<N; i++)
        {
            x[i] = xlb + ((REAL_TYPE)i)*dx;
        }

        for (int i = 0; i < N; i++)
        {
            mean_x = rho*x[i];
            for (int j = 0; j < N; j++)
            {
                if (j==0)
                {
                    x2 = x[0] + 0.5*dx;
                    Px[i*N + j] = normcdf(x2 - mean_x, sigma);
                }
                else if (j==N-1)
                {
                    x1 = x[N-1] - 0.5*dx;
                    Px[i*N + j] = 1 - normcdf(x1 - mean_x, sigma);
                }
                else
                {
                    x1 = x[j] - 0.5*dx;
                    x2 = x[j] + 0.5*dx;
                    Px[i*N + j] = normcdf(x2 - mean_x, sigma) - normcdf(x1 - mean_x, sigma);
                }
            }
        }

        #ifdef DEBUG_IDENTITY_MATRIX

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (i==j)
                {
                    Px[i*N + j] = 1.0;
                }
                else
                {
                    Px[i*N + j] = 0.0;
                }
            }
        }

        #endif

        // check Px is a valid probability matrix
        bool bValid = true;
        REAL_TYPE sum;

        for (int i = 0; i < N; i++)
        {
            sum = 0;
            for (int j = 0; j < N; j++)
            {
                if (Px[i*N+j]<0)
                {
                    mexPrintf("Tauchen: P(%d,%d)=%g<0.\n", i, j, Px[i * N + j]);
                    bValid = false;
                }
                sum += Px[i*N+j];
            }
#ifdef USE_SINGLE
            if (abs(1.0 - sum)>1e-6)
#else
            if (abs(1.0 - sum) > 1e-10)
#endif
            {
                mexPrintf("Tauchen: abs(rowsum-1)=%g.\n", abs(1.0 - sum));
                bValid = false;
            }
        }

        if (bValid)
        {
            return (EXIT_SUCCESS);
        }
        else
        {
            mexPrintf("Tauchen: did not produce valid transition probabilities.\n");
            return EXIT_FAILURE;
        }

    }
    else
    {
        mexPrintf("Tauchen: check input parameters.\n");
        return EXIT_FAILURE;
    }

}