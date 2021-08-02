// Wrappers around some of the BLAS routines that ship with MATLAB.
//#include "/usr/local/bin/MATLAB_R2009a/extern/include/blas.h"
//#include "/usr/local/bin/MATLAB_R2009a/extern/include/lapack.h"

//#include "/Applications/MATLAB_R2009b.app/extern/include/blas.h"
//#include "/Applications/MATLAB_R2009b.app/extern/include/lapack.h"

#include <cstddef>
#include "blas.h"
#include "lapack.h"


static void
daxpy_wrapper(long n, double a, double *x, long incx, double *y, long incy)
{
	daxpy_(&n, &a, x, &incx, y, &incy);
}
#undef daxpy
#define daxpy daxpy_wrapper

static void
dcopy_wrapper(long n, double *x, long incx, double *y, long incy)
{
	dcopy_ (&n, x, &incx, y, &incy);
}
#undef dcopy
#define dcopy dcopy_wrapper

static double
idamax_wrapper(long n, double* x, long incx){
	idamax_(&n, x, &incx);
}
#undef idamax
#define idamax idamax_wrapper


static double
dnrm2_wrapper(long len, double *x, long ldx){
	return dnrm2_ (&len, x, &ldx);
}
#undef dnrm2
#define dnrm2 dnrm2_wrapper

static double 
dasum_wrapper(long len, double *x, long ldx){
	return dasum_ (&len, x, &ldx);
}
#undef dasum
#define dasum dasum_wrapper


static void
dgemm_wrapper (char ta, char tb, long m, long n, long k,
	double alpha, double *a, long lda, double *b, long ldb,
	double beta, double *c, long ldc)
{
	dgemm_ (&ta, &tb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}
#undef dgemm
#define dgemm dgemm_wrapper

static void
dgemv_wrapper (char ta, long m, long n, 
	double alpha, double *a, long lda, double *b, 
	long incb, double beta, double *c, long ldc)
{
	dgemv_ (&ta, &m, &n, &alpha, a, &lda, b, &incb, &beta, c, &ldc);
}
#undef dgemv
#define dgemv dgemv_wrapper

static void
	dsyevx_wrapper(char jobz, char range, char uplo, long n, 
	double *a, long lda, double vl, double vu, long il, long iu, char *abstol, 
	double *eig)
{	
	long out_m;
	double *out_w;
	double *out_z;
	
	double *work;
	long lwork = -1;
	long *iwork = new long[5*n];
	long *ifail = new long[n];
	long info;
	long one = 1;
		
	double optimal_lwork;
	work = &optimal_lwork;
	out_w = &optimal_lwork;
	out_z = &optimal_lwork;
	
	// the following variables are supposed to be (ptrdiff_t *): 
	// &n, &lda, &il, &iu, &m, &ldz, &lwork, &iwork, &ifail, &info	
	double s = dlamch(abstol);
	dsyevx_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &s, &out_m, out_w, out_z, &one, work, &lwork, iwork, ifail, &info);
	
	lwork = optimal_lwork;
	work = new double[lwork];

//	prinf("lwork: %d", lwork);


	dsyevx_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &s, &out_m, out_w, out_z, &one, work, &lwork, iwork, ifail, &info);
	*eig = out_w[0];

	delete [] work;
	delete [] iwork;
	delete [] ifail;
//	delete [] out_w;
//	delete [] out_z;
}
#undef dsyevx
#define dsyevx dsyevx_wrapper


static void
dger_wrapper(long m, long n, double alpha, double *x, long incx, double *y, long incy, double *a, long lda)
{
	dger_ (&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
}
#undef dger
#define dger dger_wrapper


static void
dgetrf_wrapper (long m, long n, double *a, long lda, long *ipiv, long *info)
{
	dgetrf_ (&m, &n, a, &lda, ipiv, info);
}
#undef dgetrf
#define dgetrf dgetrf_wrapper

static void
dgetri_wrapper (long n, double *a, long lda, long *ipiv, long *info)
{
	double optimal_lwork;
	double *work = &optimal_lwork;

	long lwork = -1;
	dgetri_ (&n, a, &lda, ipiv, work, &lwork, info);

	lwork = optimal_lwork;
//	work = new double[lwork*lwork];
	work = new double[lwork];
	
	dgetri_ (&n, a, &lda, ipiv, work, &lwork, info);

	delete [] work;
}
#undef dgetri
#define dgetri dgetri_wrapper

static double 
ddot_wrapper(long len, double *x, long ldx, double *y, long ldy){
	return ddot_ (&len, x, &ldx, y, &ldy);
}
#undef ddot
#define ddot ddot_wrapper