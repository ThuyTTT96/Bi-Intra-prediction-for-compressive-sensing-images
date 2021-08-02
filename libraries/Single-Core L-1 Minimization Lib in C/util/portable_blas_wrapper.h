#ifndef __PORTABLE_BLAS_WRAPPER_H__

#define __PORTABLE_BLAS_WRAPPER_H__ 


// Include this file to include a standardized (cpu) blas interface.
#pragma once

#ifdef BLAS_IMPLEMENTATION_MATLAB
//#warning "Using BLAS_IMPLEMENTATION_MATLAB..."
#include "matlab_blas_wrappers.cpp"

#elif defined(BLAS_IMPLEMENTATION_SYSTEM)
//#warning "Using BLAS_IMPLEMENTATION_SYSTEM..."
#include <cblas.h>
#undef dgemm
#define dgemm dblas_sgemm
#undef dlamch
#define dlamch dblas_slamch
#undef dsyevx
#define dsyevx dblas_ssyevx

#elif defined(BLAS_IMPLEMENTATION_MKL)
//#warning "Using BLAS_IMPLEMENTATION_MKL..."
#include "mkl.h"

#elif defined(BLAS_IMPLEMENTATION_ACML)
//#warning "Using BLAS_IMPLEMENTATION_ACML..."
#include "acml_blas_wrappers.cpp"

#else
#error "Preprocessor variable BLAS_IMPLEMENTATION_* is missing!"
#endif

#endif

