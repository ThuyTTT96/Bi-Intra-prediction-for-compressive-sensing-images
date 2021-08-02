/*% MEX/C code for PaLM fast l1-minimization

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

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "portable_blas_wrapper.h"
#include "SolvePALM.h"

double max(double a, double b){
	if (a>b) return a;
	return b;
}


int sign(double a){
	if (a>=0) return 1;
	return -1;
}

#if !defined (COMPILE_MEX)
//#  include <acml.h>
//#include <stdio.h>
#include <blas.h>
#include <lapack.h>
#else
#include <mex.h>
#include <blas.h>
#include <lapack.h>

extern "C" void
mexFunction (int nl, mxArray *pl[], int nr, mxArray const *pr[])
{
	if (nr < 4){
	  mexErrMsgTxt ("[x e nIter nIter_in_total] = SolvePALM(b, A, tol, tol_int, maxIter, maxIter_alt, stop, xG)");
	}

	double *b = (double*)mxGetPr (pr[0]);
	double *A = (double*)mxGetPr (pr[1]);
	int m = mxGetM (pr[1]);
	int n = mxGetN (pr[1]);
	
	if (mxGetM (pr[0]) * mxGetN (pr[0])  != m){
		printf("mxGetM (pr[0]) * mxGetN (pr[0]) : %20.20f\n", mxGetM (pr[0]) * mxGetN (pr[0]));
		printf("A: %d x %d\n", m, n);
		mexErrMsgTxt ("SolveDALM_fast: min |x|1 + |e|1 s.t. Ax + e = b\n");
	}
	double tol = mxGetScalar (pr[2]);
	double tol_int = mxGetScalar (pr[3]);
	int maxIter = (int)mxGetScalar(pr[4]);
	int maxIter_apg = (int)mxGetScalar(pr[5]);
	double *xG = NULL;	
	int stop = 5;

	//if (nr > 6)
	  //  xG = (double *)mxGetPr(pr[7]);
	//if(nr > 5)
	//  stop = (int)mxGetScalar(pr[6]);
	
	double *x = new double[n];
	double *e = new double[m];
	int nIter;
	int nIter_in_total;
	
	SolvePALM(x, e, nIter, nIter_in_total, b, A, tol, tol_int, maxIter, maxIter_apg, m, n);
	
	if (nl > 0){
		pl[0] = mxCreateNumericMatrix (n, 1, mxDOUBLE_CLASS, mxREAL);
		memcpy (mxGetData (pl[0]), (void*)x, n*sizeof(double));
	}
	if (nl > 1){
		pl[1] = mxCreateNumericMatrix (m, 1, mxDOUBLE_CLASS, mxREAL);
		memcpy (mxGetData (pl[1]), (void*)e, m*sizeof(double));
	}
		
	if (nl > 2)	 pl[2] = mxCreateDoubleScalar (nIter);
	if (nl > 3)	 pl[3] = mxCreateDoubleScalar (nIter_in_total);
	
	delete [] x;
	delete [] e;
}

#endif

double const eps = 1e-8;

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
SolvePALM (
	double *&x, double *&e, int&  nIter, int &nIter_in_total, 
	double *b, double *A, double tol, double tol_int, int maxIter, int maxIter_apg, int m, int n){
		
		bool converged_main, converged_apg;
		
		int nIter_apg;
		int ldA = m;
		
		int maxmn = m;
		if (maxmn < n)	maxmn = n;
		
		double tau, tauInv, mu, muInv, muTauInv;
		double *lambda = new double[m];
		double *lambdaScaled = new double[m];
		//x      = new double[n];
		double *z      = new double[n];
		//e      = new double[m];
		double *x_old_main = new double[n];
		double *e_old_main = new double[m];
		double *x_old_apg  = new double[n];
		double *temp  = new double[maxmn];
		double *temp1 = new double[maxmn];
		double *temp2 = new double[maxmn];
		double t1, t2;
		
		double *G      = new double[n*n];
		double *Gx     = new double[n];
		double *Gx_old = new double[n];		
		double *Gz     = new double[n];
		double *s      = new double[n];
		
		double n1, n2, e1, e2;
		
		for(int i = 0 ; i < m; i++)  lambda[i] = 0;
		for(int i = 0 ; i < n; i++)  x[i] = 0;
		for(int i = 0 ; i < m; i++)  e[i] = b[i];
		
		dgemm ('T', 'N', n, n, m, (double)1.0, A, m, A, m, (double)0.0, G, n); 
		char ss = 'S';	
		dsyevx('N', 'I', 'U', n, G, n, 0.0 ,0.0 , n,n, &ss, &tau);
		tau = tau;
		dgemm ('T', 'N', n, n, m, (double)1.0, A, m, A, m, (double)0.0, G, n); 
		
		tauInv = 1.0 / tau;
		nIter  = 0;
		mu     = 2 * m / dasum(m, b, 1);
		converged_main = false;
		nIter_in_total = 0;
		
		
		while (!converged_main){
			
			muInv = 1.0 / mu;
			
			for(int i = 0 ; i < m; i++)  lambdaScaled[i] = muInv * lambda[i];
			nIter = nIter + 1;
			
			for(int i = 0 ; i < n; i++)  x_old_main[i] = x[i];
			for(int i = 0 ; i < m; i++)  e_old_main[i] = e[i];
			
			for(int i = 0 ; i < m; i++)  {
				temp2[i] = b[i] + lambdaScaled[i];
				temp[i]  = temp2[i];
			}
			
			dgemv('N', m, n, -1.0, A, ldA, x, 1, 1, temp, 1);
			for(int i = 0 ; i < m; i++){
				e[i] = max(fabs(temp[i]) - muInv, 0);
				if (temp[i] < 0) e[i] = -1 * e[i];
			}
			
			converged_apg = false;
			
			for(int i = 0 ; i < m; i++) temp2[i] = e[i] - temp2[i];
			dgemv('T', m, n, 1.0, A, ldA, temp2, 1, 0, temp1, 1);
			
			nIter_apg = 0;
			
			t1 = 1;
			for(int i = 0 ; i < n; i++)  z[i] = x[i];
			
			muTauInv = muInv * tauInv;
			
			dgemv('N', n, n, 1.0, G, n, x, 1, 0, Gx, 1);
			for(int i = 0 ; i < n; i++) Gz[i] = Gx[i];

			while (!converged_apg){
				nIter_apg = nIter_apg + 1;
				
				for(int i = 0 ; i < n ; i++) {
					x_old_apg[i] = x[i];
					Gx_old[i]    = Gx[i];
					temp[i]      = z[i] - tauInv * (temp1[i] + Gz[i]);
					x[i]         = (1.0*sign(temp[i])) * max(fabs(temp[i]) - muTauInv, 0);
				}				
				//Gx = G * x;
				dgemv('N', n, n, 1.0, G, n, x, 1, 0, Gx, 1);
				
				//s = tau * (z - x) + Gx - Gz;
				for(int i = 0 ; i < n; i++){
					s[i] = tau * (z[i] - x[i]) + Gx[i] - Gz[i];
				}
				
				//if norm(s) < tol_int * tau * max(1,norm(x))
		        //    converged_apg = 1;
		        //end
				if (dnrm2(n, s, 1) < tol_int * tau * max(1, dnrm2(n, x, 1))){
					converged_apg = true;
				}
				//if nIter_apg >= maxIter_apg
		        //    converged_apg = 1 ;
		        //end
				if (nIter_apg >= maxIter_apg){
					converged_apg = 1;
				}
				
				//t2 = (1+sqrt(1+4*t1*t1))/2 ;
		        //z = x + ((t1-1)/t2)*(x-x_old_apg) ;
		        //Gz = Gx + ((t1-1)/t2) * (Gx - Gx_old);
		        //t1 = t2 ;
				t2 = (1 + sqrt(1+4*t1*t1))/2;
				for(int i = 0 ; i < n; i++){
					z[i]  = x[i] + ((t1-1)/t2)*(x[i] - x_old_apg[i]);
					Gz[i] = Gx[i] + ((t1-1)/t2)*(Gx[i] - Gx_old[i]);
				}
				t1 = t2;
			}
						
			nIter_in_total += nIter_apg;
			
			dgemv('N', m, n, 1.0, A, ldA, x, 1, 0, temp,1);
			
			e1 = dnrm2(m, e_old_main, 1);
			for(int i = 0; i < m ; i++){
				lambda[i] = lambda[i] + mu * (b[i] - temp[i] - e[i]);
				e_old_main[i] -= e[i];
			}
			e2 = dnrm2(m, e_old_main, 1);
			
			n1 = dnrm2(n, x_old_main, 1);
			for(int i = 0 ; i < n; i++){
				x_old_main[i] -= x[i];
			}
			n2 = dnrm2(n, x_old_main, 1);
			
			if (n2 <= tol * n1 && e2 <= tol * e1){
				converged_main = true;
			}
			if (! converged_main && nIter >= maxIter){
				converged_main = true;
			}
			//printf("nIter:%d x1: %f  x2: %f norm(x): %f\n", nIter, x1, x2, dnrm2(n, x, 1));
			//printf("nIter:%d e1: %f  e2: %f norm(e): %f\n", nIter, e1, e2, dnrm2(n, e, 1));
		}
		
		delete [] lambda; 		
		delete [] lambdaScaled;
		delete [] x_old_main;
		delete [] z;
		delete [] e_old_main;
		delete [] x_old_apg;
		delete [] temp;
		delete [] temp1;
		delete [] temp2;
		delete [] G;
		delete [] Gx;
		delete [] Gx_old;
		delete [] Gz;
		delete [] s;
	}
