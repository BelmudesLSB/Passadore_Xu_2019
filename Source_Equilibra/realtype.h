#ifndef REALTYPE_H
#define REALTYPE_H

// NOTE: define these in the Makefile, e.g. -DUSEDOUBLE would define USE_DOUBLE
//#define USE_DOUBLE;
//#define USE_SINGLE;

#if defined(USE_SINGLE)
// #warning "Compiling single precision code."
#elif defined(USE_DOUBLE)
// #warning "Compiling double precision code."
#else
#error "Compiler option should be either -DUSE_SINGLE or -DUSE_DOUBLE"
#endif

#ifdef USE_SINGLE
typedef float REAL_TYPE;
#define CUDART_NEG_INF -CUDART_INF_F
#define POWERFUN powf
#endif

#ifdef USE_DOUBLE
typedef double REAL_TYPE;
#define CUDART_NEG_INF -CUDART_INF
#define POWERFUN pow
#endif

#endif // !FLOATINGPOINT_H
