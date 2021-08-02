// Wrappers around some of the BLAS routines that ship with MATLAB.
//#include "/usr/local/bin/MATLAB_R2009a/extern/include/blas.h"
//#include "/usr/local/bin/MATLAB_R2009a/extern/include/lapack.h"

//#include "/Applications/MATLAB_R2009b.app/extern/include/blas.h"
//#include "/Applications/MATLAB_R2009b.app/extern/include/lapack.h"

#pragma once

#include <acml.h>

#undef long
#undef int64_t
#define long int
#define int64_t int

//#include <stdio.h>

static void
sgemm_wrapper (char ta, char tb, long m, long n, long k,
	float alpha, float *a, long lda, float *b, long ldb,
	float beta, float *c, long ldc)
{
  sgemm(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc );
}
#undef sgemm
#define sgemm sgemm_wrapper

static void
	ssyevx_wrapper(char jobz, char range, char uplo, int64_t n, 
	float *a, int64_t lda, float vl, float vu, int64_t il, int64_t iu, char *abstol, 
	float *eig)
{	
	int64_t out_m=1;
	float *out_w;
	float *out_z;
	
	float *work;
	//	int64_t lwork = -1;
	int64_t *iwork = new int64_t[5*n];
	int64_t *ifail = new int64_t[n];
	int64_t info;
	int64_t one = 1;
		
	float optimal_lwork;
	work = &optimal_lwork;
	out_w = &optimal_lwork;
	out_z = &optimal_lwork;
	
	//printf("here");

	// the following variables are supposed to be (ptrdiff_t *): 
	// &n, &lda, &il, &iu, &m, &ldz, &lwork, &iwork, &ifail, &info	
	ssyevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, slamch(*abstol), &out_m, out_w, out_z, one,ifail, &info);
	
	//	lwork = optimal_lwork;
	//	work = new float[lwork];

//	prinf("lwork: %d", lwork);


	ssyevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, slamch(*abstol), &out_m, out_w, out_z, one, ifail, &info);
	*eig = out_w[0];

	delete [] work;
	delete [] iwork;
	delete [] ifail;
//	delete [] out_w;
//	delete [] out_z;
}
#undef ssyevx
#define ssyevx ssyevx_wrapper

static void
sgetrf_wrapper (int64_t m, int64_t n, float *a, int64_t lda, int64_t *ipiv, int64_t *info)
{
	sgetrf_ (&m, &n, a, &lda, ipiv, info);
}
#undef sgetrf
#define sgetrf sgetrf_wrapper

static void
sgetri_wrapper (int64_t n, float *a, int64_t lda, int64_t *ipiv, int64_t *info)
{
	float optimal_lwork;
	float *work = &optimal_lwork;

	long lwork = -1;
	sgetri_ (&n, a, &lda, ipiv, work, &lwork, info);

	lwork = optimal_lwork;
//	work = new float[lwork*lwork];
	work = new float[lwork];
	
	sgetri_ (&n, a, &lda, ipiv, work, &lwork, info);

	delete [] work;
}
#undef sgetri
#define sgetri sgetri_wrapper
