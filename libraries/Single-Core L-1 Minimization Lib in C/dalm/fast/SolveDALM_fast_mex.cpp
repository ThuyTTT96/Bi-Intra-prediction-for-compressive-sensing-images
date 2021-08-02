/*% MEX/C code for DALM fast l1-minimization

% Copyright ©2010. The Regents of the University of California (Regents).
% All Rights Reserved. Contact The Office of Technology Licensing,
% UC Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA 94720-1620,
% (510) 643-7201, for commercial licensing opportunities.

% Created by Victor Shia, Mark Murphy, Allen Y. Yang, Department of EECS, University of California,
% Berkeley.

% IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
% SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
% ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
% REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
% PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY,
% PROVIDED HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO
% PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
*/

#include <string.h>
#include <math.h>
#include "portable_blas_wrapper.h"
#include "SolveDALM_fast.h"

int max(int a, int b){
	if (a>b) return a;
	return b;
}

#if !defined (COMPILE_MEX)
#  include <stdio.h>
#else
#  include <mex.h>

#undef printf
#define printf mexPrintf

extern "C" void
mexFunction (int nl, mxArray *pl[], int nr, mxArray const *pr[])
{
	if (nr < 4){
	  mexErrMsgTxt ("[x nIter] = solveDALM_fast(A, b, nu, tol, stop, xG)");
	}

	double *b = mxGetPr (pr[1]);
	double *A = mxGetPr (pr[0]);
	int m = mxGetM (pr[0]);
	int n = mxGetN (pr[0]);
	

	if (mxGetM (pr[1]) * mxGetN (pr[1])  != m){
		mexErrMsgTxt ("solveDALM_fast: min |x|1 + |e|1 s.t. Ax + e = b\n");
	}

	double nu = mxGetScalar (pr[2]);
	double tol = mxGetScalar (pr[3]);
	int maxIter = (int)mxGetScalar(pr[4]);
	double *xG = NULL;	
	int stop = 5;

	if (nr > 6)
	  xG = mxGetPr(pr[6]);
	if(nr > 5)
	  stop = (int)mxGetScalar(pr[5]);
	
	double *x = new double[n];
	int nIter;

	SolveDALM_fast(x, nIter, b, A, nu, tol, maxIter, m, n, stop, xG);
	//	norm_x_e_dual (x, e, nIter, b, A, nu, tol, m, n);

	if (nl > 0){
		pl[0] = mxCreateNumericMatrix (n, 1, mxDOUBLE_CLASS, mxREAL);
		memcpy (mxGetData (pl[0]), (void*)x, n*sizeof(double));
	}
	delete [] x;

	if (nl > 1){
		pl[1] = mxCreateDoubleScalar (nIter);
	}
	
	//pl[0] = mxCreateDoubleScalar(10);
}

#endif

double const eps = 1e-15;

enum stoppingCriteria{
  STOPPING_GROUND_TRUTH   = -1,
  STOPPING_DUALITY_GAP    = 1,
  STOPPING_SPARSE_SUPPORT = 2,
  STOPPING_OBJECTIVE_VALUE = 3,
  STOPPING_SUBGRADIENT    = 4,
  STOPPING_INCREMENTS     = 5,
  STOPPING_DEFAULT        = STOPPING_INCREMENTS
};

void
SolveDALM_fast (
	double *&x, int&  nIter,
	double *b, double *A, double lambda, double tol, int maxIter, int m, int n, int stoppingCriterion, double *xG)
{
	int ldA = m;

	enum stoppingCriteria stop;
	switch (stoppingCriterion){
		case -1:
		stop = STOPPING_GROUND_TRUTH;
		break;
		case 1:
		stop = STOPPING_DUALITY_GAP;
		break;
		case 2:
		stop = STOPPING_SPARSE_SUPPORT;
		break;
		case 3:
		stop = STOPPING_OBJECTIVE_VALUE;
		break;
		case 4:
		stop = STOPPING_SUBGRADIENT;
		break;
		case 5:
		stop = STOPPING_INCREMENTS;
		break;
	}

	bool verbose = false;

// beta = norm(b,1) / m;
	double beta = 0;
	for (long i = 0; i < m; i++){
		beta += fabs (b[i]);
	}
	beta = beta / m;

	double betaInv = 1 / beta;

	nIter = 0;

// y = zeros(m,1)
	double *y = new double[m];
	for (long i = 0; i < m; i++)
		y[i] = 0;

// x = zeros(n,1)
	//x = new double[n];
	for (long i = 0; i < n; i++)
		x[i] = 0;

// z = zeros (m+n,1);
	double *z = new double[n]; 
	for (long i = 0; i < n; i++)
		z[i] = 0;

	bool converged_main = false;

// temp = A' * y;
	double *temp = new double [max(m,n)];
	dgemv ('T', m, n, 1.0, A, ldA, y, 1, 0.0, temp, 1);
 
	double *x_old  = new double[n];
	double *temp1 = new double[max(m,n)];
	double *tmp = new double[max(m,n)];
	double *g = new double[m];

// f = norm(x,1);  x is 0 at this point
	double f = 0;
	double prev_f = 0;
	double total = 0;
	double nxo, nx, dx;
	
	double *Ag = new double[n];
	
	double dg, dAg, alpha;
	
	do {
		nIter = nIter + 1;
		if(verbose) printf("==== [%d] ====\n", nIter);

	// x_old = x
		for (long i = 0; i < n; i++)
			x_old[i] = x[i];

	// % update z
	// temp1 = temp + x * betaInv
	// z = sign(temp1) .* min(1, abs(temp1));
		for (long i = 0; i < n; i++){
			temp1[i] = temp[i] + x[i] * betaInv;
		}
		for (long i = 0; i < n; i++){
			z[i] = (temp1[i] > 0 ? 1 : -1)
			       * ((fabs(temp1[i]) > 1) ? 1 : fabs(temp1[i]));
		}

	//    %compute A' * y    
	//    g = lambda * y - b + A * (beta * (temp - z) + x);
		for(long i = 0 ; i < n; i++){
			tmp[i] = beta * (temp[i] - z[i]) + x[i];
		}
		for(long i = 0; i < m; i++){
			g[i] = lambda * y[i] - b[i];
		}
		dgemv('N', m, n, 1.0, A, ldA, tmp, 1, 1, g, 1);
	
	//    %alpha = g' * g / (g' * G * g);
	//    Ag = A' * g;
		dgemv('T', m, n, 1.0, A, ldA, g, 1, 0.0, Ag, 1);
	
	//    alpha = g' * g / (lambda * g' * g + beta * Ag' * Ag);
		dg = 0;
		dAg = 0;
		for(long i = 0 ; i < n; i++){
			dAg += Ag[i] * Ag[i];
		}
		for(long i = 0 ; i < m; i++){
			dg += g[i] * g[i];
		}
		alpha = dg / (lambda * dg + beta * dAg);
	
	//    y = y - alpha * g;
		for(long i = 0 ; i < m; i++){
			y[i] = y[i] - alpha * g[i];
		}
	
	//    temp = A' * y;
		dgemv ('T', m, n, 1.0, A, ldA, y, 1, 0.0, temp, 1);

	// % update x
	// x = x - beta * (z - temp);
		for (long i = 0; i < n; i++)
			x[i] = x[i] - beta * (z[i] - temp[i]);

		switch (stop){
		case STOPPING_GROUND_TRUTH:
		  total = 0;
		  for(int i = 0 ; i < n; i++){
		    total += (xG[i] - x[i])*(xG[i] - x[i]);
		  }
		  if (total < tol * tol)
		    converged_main = true;
		  break;
		case STOPPING_SUBGRADIENT:
		  printf("Duality gap is not a valid stopping criterion for ALM.");
		  break;
		case STOPPING_SPARSE_SUPPORT:
		  printf("DALM does not have a support set.");
		  break;
		case STOPPING_OBJECTIVE_VALUE:
		  prev_f = f;
		  f = 0;
		  for(int i = 0 ; i < n; i++){
		    f += fabs(x[i]);
		  }
		  if (fabs(f-prev_f)/prev_f <= tol){
		    converged_main = true;
		  }
		  break;
		case STOPPING_DUALITY_GAP:
		  printf("Duality gap is not a valid stopping criterion for ALM.");
		  break;
		case STOPPING_INCREMENTS:
		  // if norm(x_old - x) < tol * norm(x_old)
		  //     converged_main = true;
		  
		  nxo = 0;
		  for (int i = 0; i < n; i++)
		    nxo = nxo + x_old[i]*x_old[i];
		  nxo = sqrt (nxo);
		  
		  nx = 0;
		  for (int i = 0; i < n; i++)
		    nx = nx + x[i]*x[i];
		  nx = sqrt (nx);
		  
		  dx = 0;
		  for (int i = 0; i < n; i++)
		    dx = dx + (x_old[i]-x[i])*(x_old[i]-x[i]);
		  dx = sqrt(dx);
		  
		  if (dx < tol*nxo)
		    converged_main = true;
		  
		  if (verbose){
		    printf("  ||x|| = %f\n", nx);
		  }
		  
		  if (verbose){
		    if (nIter > 1){
		      printf ("  ||dx|| = %f (= %f * ||x_old||)\n",
				 dx, dx/(nxo+eps));
		    } else {
		      printf ("  ||dx|| = %f\n", dx);
		    }
		  }
		  break;
		default:
		  printf("Undefined stopping criterion.");
		  break;
		}

 		if (nIter >= maxIter){
			if (verbose)
				printf ("Maximum Iterations Reached\n");
			converged_main = true;
		}
	}
	while (!converged_main);

	if (verbose) printf("==== CONVERGED ==== \n", nIter);

	delete [] tmp;
	delete [] g;
	delete [] Ag;
	delete [] y;
	delete [] z;
	delete [] x_old;
	delete [] temp;
	delete [] temp1;
}
