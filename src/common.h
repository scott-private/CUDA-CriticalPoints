/*
============================================================================
Author      : Scott Fu
Date        : 17/06/2020
Copyright   : scottfu@foxmail.com
File Name   : common.h
============================================================================
*/

#ifndef COMMON_H
#define COMMON_H

#include <cuda_runtime.h>
#include <curand_kernel.h>

#define EXACT                           // Turns on/off exact arithmetics
#define SINGLE_PRECISION

typedef unsigned int uint; 

#ifdef SINGLE_PRECISION
#define REAL float
#else
#define REAL double
#endif

#define EPS    1e-9

#define MINREAL -1E10
#define MAXREAL 1E10
#define MAXINT 0x7FFFFFFF
#define MININT -0x7FFFFFFF
#define PI 3.141592654



#define checkCudaErrors(err) __checkCudaErrors(err, __FILE__, __LINE__)

// These are the inline versions for all of the SDK helper functions
inline void __checkCudaErrors(cudaError_t err, const char *file, const int line) {
  if (cudaSuccess != err) {
    fprintf(stderr,
            "checkCudaErrors() Driver API error = %04d \"%s\" from file <%s>, "
            "line %i.\n",
            err, cudaGetErrorString(err), file, line);
    exit(EXIT_FAILURE);
  }
}


#endif
