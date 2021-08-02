// % MEX/C code for Homotopy l1-minimization
// 
// % Copyright ©2010. The Regents of the University of California (Regents).
// % All Rights Reserved. Contact The Office of Technology Licensing,
// % UC Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA 94720-1620,
// % (510) 643-7201, for commercial licensing opportunities.
// 
// % Created by Victor Shia, Mark Murphy, Allen Y. Yang, Department of EECS, University of California,
// % Berkeley.
// 
// % IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
// % SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
// % ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
// % REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// % REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED
// % TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// % PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY,
// % PROVIDED HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO
// % PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.


#include <string.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <cstdlib>
#include "portable_blas_wrapper.h"
#include "SolveHomotopy.h"

int sign(double a){
	if (a > 0.0){
		return 1;
	}else{
		return -1;
	}
}

static int
min_wrapper(int a, int b){
	return a <= b ? a : b;
}
#undef min
#define min min_wrapper


#if !defined (COMPILE_MEX)
#include <stdio.h>
#else
#include <mex.h>
#undef printf
#define printf mexPrintf
//#undef long
//#define long int

extern "C" void
mexFunction (int nl, mxArray *pl[], int nr, mxArray const *pr[])
{
	if (nr < 4){
	  mexErrMsgTxt ("[x nIter] = SolveHomotopy(A, b, nu, tol, stop, xG)");
	}

	double *b = mxGetPr (pr[1]);
	double *A = mxGetPr (pr[0]);
	int m = mxGetM (pr[0]);
	int n = mxGetN (pr[0]);
	int i = 0;

	double nu = mxGetScalar (pr[2]);
	double tol = mxGetScalar (pr[3]);
	int maxIter = (int)mxGetScalar(pr[4]);
//	printf("tolerance: %f\n", tol);
	double *xG;	int stop;

	if (mxGetM (pr[1]) * mxGetN (pr[1])  != m){
		mexErrMsgTxt ("SolveHomotopy: min |x|1 + |e|1 s.t. Ax + e = b\n");
	}
	
	if (nr < 6)
	  xG = NULL;
	else
	  xG = mxGetPr(pr[6]);
	if(nr < 5)
	  stop = 3;
	else
	  stop = (int)mxGetScalar(pr[5]);
	
	int nIter;
	
	double *x;
	SolveHomotopy(x, nIter, b, A, tol, nu, maxIter, m, n);
	
	//nl = 1;
	//pl[0] = mxCreateDoubleScalar(0);
	
	if (nl > 0){
		//pl[0] = mxCreateDoubleScalar(0);
		pl[0] = mxCreateNumericMatrix (n, 1, mxDOUBLE_CLASS, mxREAL);
		memcpy (mxGetData (pl[0]), (void*)x, n*sizeof(double));
	}
	
	delete [] x;

	if (nl > 1){
		pl[1] = mxCreateDoubleScalar (nIter);
	}
	
}

#endif

double *getMatrixFromIndices(double *A, int m, int n, int *index, int len){
	double *rtn = new double[m * len];
	for(int i = 0 ; i < len; i++){
		for(int j = 0; j < m; j++){
			rtn[i * m + j] = A[index[i]*m + j];
		}
	}
	return rtn;
}

// global N gamma_x z_x  xk_temp del_x_vec pk_temp dk epsilon isNonnegative
int N;
int *gamma_x;
int gamma_x_size;
int *z_x;
double *xk_temp;
double *del_x_vec;
double *pk_temp;
double *dk;
double epsilon;
//Real eps = D_MACHEPS;
double eps = 2.2e-16;
//double tolerance = 0.0001;
//double lambda = 1e-4;

void update_primal(int &i_delta, double &delta,int &out_x);

double* SolveHomotopy(double *&x, int &nIter, double *y, double *A, double tol, double lambda, int maxIter, int m, int n){
	int lda = m;
	N = n;
//	lambda = 1e-3;
//	maxIter = 100;
//	isNonnegative = true;
//	verbose = false;
//	xk_1 = [];
	bool verbose = false;
	double *xk_1 = new double[n];
	double *x_k  = new double[n];
	memset(xk_1, 0, n * sizeof(double));
 	memset(x_k, 0, n * sizeof(double));
	
//	z_x = zeros(N,1);
//	gamma_x = [];       % Primal support
	z_x = new int[n];
	memset(z_x, 0, n * sizeof(int));
	
	int *z_xk = new int[n];
	memset(z_xk, 0, n * sizeof(int));
	
	gamma_x = new int[min(maxIter, n)];
	memset(gamma_x, 0, min(maxIter, n)*sizeof(int));
	gamma_x_size = 0;
	
	int *gamma_xk = new int[min(maxIter, n)];
	memset(gamma_xk, 0, min(maxIter, n)*sizeof(int));
	int gamma_xk_size = 0;
	
//	Primal_constrk = -A'*y;
	double *Primal_constrk = new double[n];
//	memset(Primal_constrk, 0, n*sizeof(double));
	pk_temp = new double [n];
	dgemv('T', m, n, -1, A, m, y, 1, 0, Primal_constrk, 1);
	
//    [c i] = max(abs(Primal_constrk));	
	double c = idamax(n, Primal_constrk, 1);
	c = idamax(n, Primal_constrk, 1); // for some reason, I need to call idamax twice or else it doesn't work.
	int i = -1;
	for(int k = 0; k < n; k++){
		if (Primal_constrk[k] == c || Primal_constrk[k] == -1 * c) i = k;
	}
//	int i = 1000;
	
	epsilon = c;
	
//	nz_x = zeros(N,1);
	int *n_zx = new int[n];
	
//	xk_1 = zeros(N,1);
//  gamma_xk = i;
	gamma_xk[gamma_xk_size++] = i;
	
//	f = epsilon*norm(xk_1,1) + 1/2*norm(y-A*xk_1)^2;
//	z_x(gamma_xk) = -sign(Primal_constrk(gamma_xk));
	double *tempD = new double[m];
	dcopy(m, y, 1, tempD, 1);
	dgemv('N', m, n, -1, A, lda, xk_1, 1, 1, tempD, 1);
	double prev_f;
	double f = epsilon * dasum(n, xk_1, 1) + 0.5 * pow(dnrm2(m, tempD, 1), 2);
	delete [] tempD;
	
	for(int k = 0; k < gamma_xk_size; k++){
		z_x[gamma_xk[k]] = -1 * sign(Primal_constrk[gamma_xk[k]]);
	}
	
//	z_xk = z_x;
//	dcopy(n, z_x, 1, z_xk, 1);
	memcpy(z_xk, z_x, n*sizeof(int));
	
//	iter = 0;
//	out_x = [];
//	old_delta = 0;
//	count_delta_stop = 0;
	int iter = 0;
	int out_x;
	double old_delta = 0;
	int count_delta_stop = 0;

//	AtgxAgx = A(:,gamma_xk)'*A(:,gamma_xk);
//	iAtgxAgx = inv(A(:,gamma_xk)'*A(:,gamma_xk));
//  double *mTemp1 = getMatrixFromIndices(A, 0, m, gamma_xk, 1); // NEED TO WRITE
	double *AtgxAgx = new double[maxIter * maxIter];
	int ldAtgxAgx = maxIter;
	
	AtgxAgx[0] = ddot(m, (A + lda * gamma_xk[0]), 1, (A + lda * gamma_xk[0]), 1);
		
	double *iAtgxAgx = new double[maxIter * maxIter];
	int AtgxAgx_dim = 1;
	iAtgxAgx[0] = 1.0 / AtgxAgx[0];
	
	double *del_x;
	del_x_vec = new double[n];
	memset(del_x_vec, 0, n*sizeof(double));
	
	double *Asupported;
	double *Agdelx = new double[m];;
	dk = new double[n];
	
	double minEps = epsilon < eps * 2 ? epsilon : eps * 2;//min(epsilon, eps*2);
	bool keep_going;
	double *tempD2;
	int len_gamma;
	int outx_index;
	int rowi, colj;
	int temp_int;
	double *Q12, *Q21, Q22;
	double *Q12Q21_Q22;
	double *AtgxAnx;
	double *iA11A12, *A21iA11;
	double S;
	double *Q11_right;
	int i_delta;
	double delta;
	double epsilon_old;
	double tempFloat;
	xk_temp = new double[n];

	while (iter < maxIter){
//		iter = iter+1;	    
		iter++;
		
//		gamma_x = gamma_xk;
//	    z_x = z_xk;
//	    x_k = xk_1;
//		dcopy(gamma_xk_size, gamma_xk, 1, gamma_x, 1);
//		dcopy(n, z_xk, 1, z_x, 1);
		memcpy(gamma_x, gamma_xk, gamma_xk_size * sizeof(int));
		gamma_x_size = gamma_xk_size;
		memcpy(z_x, z_xk, sizeof(int)*n);
		dcopy(n, xk_1, 1, x_k, 1);
		
//		%%%%%%%%%%%%%%%%%%%%%
//	    %%%% update on x %%%%
//	    %%%%%%%%%%%%%%%%%%%%%
		
//	    % Update direction
//	    del_x = iAtgxAgx*z_x(gamma_x);
//	    del_x_vec = zeros(N,1);
//	    del_x_vec(gamma_x) = del_x;
		del_x = new double[gamma_x_size];
		tempD = new double[gamma_x_size];
		for(int k = 0; k < gamma_x_size; k++){
			tempD[k] = z_x[gamma_x[k]];
		}
		dgemv('N', AtgxAgx_dim, AtgxAgx_dim, 1, iAtgxAgx, ldAtgxAgx, tempD, 1, 0, del_x, 1);
		memset(del_x_vec, 0, N*sizeof(double));
		delete [] tempD;
		
		for(int k = 0; k < gamma_x_size; k++){
			del_x_vec[gamma_x[k]] = del_x[k];
		}
		
//	    %dk = A'*(A*del_x_vec);
//	    Asupported = A(:,gamma_x);
//	    Agdelx = Asupported*del_x;
//	    dk = A'*Agdelx;
		Asupported = getMatrixFromIndices(A, m, n, gamma_x, gamma_x_size);		
		dgemv('N', m, gamma_x_size, 1, Asupported, m, del_x, 1, 0, Agdelx, 1);
		dgemv('T', m, n, 1, A, lda, Agdelx, 1, 0, dk, 1);
		
//		%%% CONTROL THE MACHINE PRECISION ERROR AT EVERY OPERATION: LIKE BELOW. 
//	    pk_temp = Primal_constrk;
//	    gammaL_temp = find(abs(abs(Primal_constrk)-epsilon)<min(epsilon,2*eps));
//	    pk_temp(gammaL_temp) = sign(Primal_constrk(gammaL_temp))*epsilon;
		dcopy(n, Primal_constrk, 1, pk_temp, 1);
		minEps = epsilon < eps * 2 ? epsilon : eps * 2;
		for(int k = 0 ; k < n; k++){
			if (fabs(fabs(Primal_constrk[k])-epsilon) < minEps)
				pk_temp[k] = sign(Primal_constrk[k]) * epsilon;
		}
		
//		xk_temp = x_k;
//	    xk_temp(abs(x_k)<2*eps) = 0;
		dcopy(n, x_k, 1, xk_temp, 1);
		for(int k = 0; k < n; k++){
			if (fabs(x_k[k]) < 2*eps){
				xk_temp[k] = 0;
			}
		}
		
//		[i_delta, delta, out_x] = update_primal(out_x);
		update_primal(i_delta, delta, out_x);
		
/*		
		if old_delta < 4*eps && delta < 4*eps
	        count_delta_stop = count_delta_stop + 1;

	        if count_delta_stop >= 500
	            if verbose
	                disp('stuck in some corner');
	            end
	            break;
	        end
	    else
	        count_delta_stop = 0;
	    end
	    old_delta = delta;
*/		    
		if(old_delta < 4*eps && delta < 4*eps){
			count_delta_stop += 1;
			if (count_delta_stop >= 500){
				if (verbose)  printf("stuck in some corner");
				break;
			}
		}else{
			count_delta_stop = 0;
		}
		old_delta = delta;
		
//	    xk_1 = x_k+delta*del_x_vec;
//	    Primal_constrk = Primal_constrk+delta*dk;
//	    epsilon_old = epsilon;
//	    epsilon = epsilon-delta;	
		dcopy(n, x_k, 1, xk_1, 1);
		daxpy(n, delta, del_x_vec, 1, xk_1, 1);
		daxpy(n, delta, dk, 1, Primal_constrk, 1);
		epsilon_old = epsilon;
		epsilon     = epsilon - delta;
		
//		if epsilon <= lambda;
//	        xk_1 = x_k + (epsilon_old-lambda)*del_x_vec;
//	        break;
//	    end
	
		if (epsilon <= lambda || delta == -1){
			dcopy(n, x_k, 1, xk_1, 1);
			daxpy(n, (epsilon_old - lambda), del_x_vec, 1, xk_1, 1);
			break;
		}
		
//		keep_going = true;
		keep_going = true;
		
/*
		if delta~=0
	        switch stoppingCriterion
	            case STOPPING_GROUND_TRUTH
	                keep_going = norm(xk_1-xG)>tolerance;
	            case STOPPING_SPARSE_SUPPORT
	                nz_x_prev = nz_x;
	                nz_x = (abs(xk_1)>eps*10);
	                num_nz_x = sum(nz_x(:));
	                num_changes_active = (sum(nz_x(:)~=nz_x_prev(:)));
	                if num_nz_x >= 1
	                    criterionActiveSet = num_changes_active / num_nz_x;
	                    keep_going = (criterionActiveSet > tolerance);
	                end
	            case STOPPING_DUALITY_GAP
	                error('Duality gap is not a valid stopping criterion for Homotopy.');
	            case STOPPING_OBJECTIVE_VALUE
	                % continue if not yeat reached target value tolA
	                prev_f = f;
	                f = lambda*norm(xk_1,1) + 1/2*norm(y-Asupported*xk_1(gamma_x))^2;
	                keep_going = (abs((prev_f-f)/prev_f) > tolerance);
	            case STOPPING_SUBGRADIENT
	                keep_going = norm(delta*del_x_vec)>tolerance;
	            otherwise,
	                error('Undefined stopping criterion');
	        end % end of the stopping criteria switch
	    end
*/
		if (delta != 0){
			prev_f = f;
			tempD = new double[gamma_x_size];
			for(int k = 0; k < gamma_x_size; k++){
				tempD[k] = xk_1[gamma_x[k]];
			}
			tempD2 = new double[m];
			dcopy(m, y, 1, tempD2, 1);
			dgemv('N', m, gamma_x_size, -1, Asupported, m, tempD, 1, 1, tempD2, 1);
			
			f = lambda * dasum(n, xk_1, 1) + 0.5 * pow(dnrm2(m, tempD2, 1), 2);
			keep_going = fabs((prev_f - f)/prev_f) > tol;
			delete [] tempD2;
			delete [] tempD;
		}
		
//	    if keep_going && norm(xk_1 - x_k)<100*eps
//	        if verbose
//	            disp('The iteration is stuck.');
//	        end
//	        keep_going = false;
//	    end
		tempD = new double[n];
		dcopy(n, x_k, 1, tempD, 1);
		daxpy(n, -1, xk_1, 1, tempD, 1);
		if (keep_going && dnrm2(n, tempD, 1) < 100*eps){
			keep_going = false;
		}
		delete [] tempD;
		
//	    if ~keep_going
//	        break;
//	    end
		if (keep_going == false){
			break;
		}
		
		if (out_x != -1){
//			% If an element is removed from gamma_x
//	        len_gamma = length(gamma_x);
//	        outx_index = find(gamma_x==out_x(1));
//	        gamma_x(outx_index) = gamma_x(len_gamma);
//	        gamma_x(len_gamma) = out_x(1);
//	        gamma_x = gamma_x(1:len_gamma-1);
//	        gamma_xk = gamma_x;
			len_gamma = gamma_x_size;
			outx_index = -1;
			for(int k = 0 ; k < gamma_x_size; k++){
				if(gamma_x[k] == out_x){
					outx_index = k;
				}
			}
			gamma_x[outx_index] = gamma_x[len_gamma-1];
			gamma_x[len_gamma-1] = out_x;
			gamma_x_size--;
//			dcopy(gamma_x_size, gamma_x, 1, gamma_xk, 1);
			memcpy(gamma_xk, gamma_x, gamma_x_size * sizeof(int));
			gamma_xk_size = gamma_x_size;

/*			
			rowi = outx_index; % ith row of A is swapped with last row (out_x)
	        colj = outx_index; % jth column of A is swapped with last column (out_lambda)
	        AtgxAgx_ij = AtgxAgx;
	        temp_row = AtgxAgx_ij(rowi,:);
	        AtgxAgx_ij(rowi,:) = AtgxAgx_ij(len_gamma,:);
	        AtgxAgx_ij(len_gamma,:) = temp_row;
	        temp_col = AtgxAgx_ij(:,colj);
	        AtgxAgx_ij(:,colj) = AtgxAgx_ij(:,len_gamma);
	        AtgxAgx_ij(:,len_gamma) = temp_col;
*/
			rowi = outx_index;
			colj = outx_index;
					
			for(int k = 0 ; k < AtgxAgx_dim; k++){
				tempFloat = AtgxAgx[rowi + k * ldAtgxAgx];
				AtgxAgx[rowi + k * ldAtgxAgx] = AtgxAgx[len_gamma-1 + k*ldAtgxAgx];
				AtgxAgx[len_gamma-1 + k*ldAtgxAgx] = tempFloat;
			}
			for(int k = 0 ; k < AtgxAgx_dim; k++){
				tempFloat = AtgxAgx[colj * ldAtgxAgx + k];
				AtgxAgx[colj * ldAtgxAgx + k] = AtgxAgx[(len_gamma - 1) * ldAtgxAgx + k];
				AtgxAgx[(len_gamma - 1) * ldAtgxAgx + k] = tempFloat;
			}
			
/*			
	        iAtgxAgx_ij = iAtgxAgx;
	        temp_row = iAtgxAgx_ij(colj,:);
	        iAtgxAgx_ij(colj,:) = iAtgxAgx_ij(len_gamma,:);
	        iAtgxAgx_ij(len_gamma,:) = temp_row;
	        temp_col = iAtgxAgx_ij(:,rowi);
	        iAtgxAgx_ij(:,rowi) = iAtgxAgx_ij(:,len_gamma);
	        iAtgxAgx_ij(:,len_gamma) = temp_col;
*/			
			for(int k = 0 ; k < AtgxAgx_dim; k++){
				tempFloat = iAtgxAgx[colj + k * ldAtgxAgx];
				iAtgxAgx[colj + k * ldAtgxAgx] = iAtgxAgx[len_gamma-1 + k*ldAtgxAgx];
				iAtgxAgx[len_gamma-1 + k*ldAtgxAgx] = tempFloat;
			}
			for(int k = 0 ; k < AtgxAgx_dim; k++){
				tempFloat = iAtgxAgx[rowi * ldAtgxAgx + k];
				iAtgxAgx[rowi * ldAtgxAgx + k] = iAtgxAgx[(len_gamma - 1) * ldAtgxAgx + k];
				iAtgxAgx[(len_gamma - 1) * ldAtgxAgx + k] = tempFloat;
			}
//			AtgxAgx = AtgxAgx_ij(1:len_gamma-1,1:len_gamma-1);
			AtgxAgx_dim--;
			
//			n = size(AtgxAgx_ij,1);
			temp_int = AtgxAgx_dim+1;	        
			
//	        Q11 = iAtgxAgx_ij(1:n-1,1:n-1);
//	        Q12 = iAtgxAgx_ij(1:n-1,n);
//	        Q21 = iAtgxAgx_ij(n,1:n-1);
//	        Q22 = iAtgxAgx_ij(n,n);
//	        Q12Q21_Q22 = Q12*(Q21/Q22);
//	        iAtgxAgx = Q11 - Q12Q21_Q22;
			Q12 = new double[temp_int-1];
			Q21 = new double[temp_int-1];
			for(int k = 0 ; k < temp_int-1; k++){
				Q12[k] = iAtgxAgx[(temp_int-1)*ldAtgxAgx + k];
				Q21[k] = iAtgxAgx[(temp_int-1) + k*ldAtgxAgx];
			}
			Q22 = iAtgxAgx[(temp_int-1)*ldAtgxAgx + (temp_int-1)];
			//Q12Q21_Q22 = new double[(temp_int-1)*(temp_int-1)];
			//memset(Q12Q21_Q22, 0, (temp_int-1)*(temp_int-1));
			dger(temp_int-1, temp_int-1, -1/Q22, Q12, 1, Q21, 1, iAtgxAgx, ldAtgxAgx);
			
			xk_1[out_x] = 0;
			
			delete [] Q12;
			delete [] Q21;
			//delete [] Q12Q21_Q22;
			
		}else{
//			gamma_xk = [gamma_x; i_delta];
//	        new_x = i_delta;
//			dcopy(gamma_x_size, gamma_x, 1, gamma_xk, 1);
			memcpy(gamma_xk, gamma_x, gamma_x_size*sizeof(int));
			gamma_xk_size = gamma_x_size;
			gamma_xk[gamma_xk_size] = i_delta;
			gamma_xk_size++;
			
//			AtgxAnx = A(:,gamma_x)'*A(:,new_x);
			tempD = new double[m];
			for(int k = 0 ; k < m ; k++){
				tempD[k] = A[k + i_delta*lda];
			}
			tempD2 = new double[gamma_x_size * m];
			for(int k = 0 ; k < gamma_x_size; k++){
				for(int l = 0 ; l < m ; l++){
					tempD2[l + k * m] = A[l + gamma_x[k] * m];
				}
			}
			AtgxAnx = new double[gamma_x_size];
			dgemv('T', m, gamma_x_size, 1, tempD2, m, tempD, 1, 0, AtgxAnx, 1);
			delete [] tempD2;
			
//			AtgxAgx_mod = [AtgxAgx AtgxAnx; AtgxAnx' A(:,new_x)'*A(:,i_delta)];			
//			AtgxAgx = AtgxAgx_mod;
			AtgxAgx_dim++;
			for(int k = 0 ; k < gamma_x_size; k++){
				AtgxAgx[(gamma_x_size) * ldAtgxAgx + k] = AtgxAnx[k];
				AtgxAgx[gamma_x_size + k * ldAtgxAgx]   = AtgxAnx[k];
			}
			AtgxAgx[gamma_x_size * ldAtgxAgx + gamma_x_size] = ddot(m, tempD, 1, tempD, 1);
			delete [] tempD;
//			%iAtgxAgx = update_inverse(AtgxAgx, iAtgxAgx,1);
//	        n = size(AtgxAgx,1);
//			iA11 = iAtgxAgx;
//	        iA11A12 = iA11*AtgxAgx(1:n-1,n);
//	        A21iA11 = AtgxAgx(n,1:n-1)*iA11;
//	        S = AtgxAgx(n,n)-AtgxAgx(n,1:n-1)*iA11A12;
//	        Q11_right = iA11A12*(A21iA11/S);
			tempD = new double[AtgxAgx_dim-1];
			tempD2 = new double[AtgxAgx_dim-1];
			for(int k = 0 ; k < AtgxAgx_dim-1; k++){
				tempD[k]  = AtgxAgx[(AtgxAgx_dim-1)*ldAtgxAgx + k];
				tempD2[k] = AtgxAgx[k * ldAtgxAgx + AtgxAgx_dim - 1];
			}
			iA11A12 = new double[AtgxAgx_dim-1];
			A21iA11 = new double[AtgxAgx_dim-1];
			dgemv('N', AtgxAgx_dim-1, AtgxAgx_dim-1, 1, iAtgxAgx, ldAtgxAgx, tempD, 1, 0, iA11A12, 1);
			dgemv('T', AtgxAgx_dim-1, AtgxAgx_dim-1, 1, iAtgxAgx, ldAtgxAgx, tempD2, 1, 0, A21iA11, 1);
			S = AtgxAgx[(AtgxAgx_dim-1)*ldAtgxAgx+(AtgxAgx_dim-1)] - ddot(AtgxAgx_dim-1, tempD2, 1, iA11A12, 1);
			Q11_right = new double[(AtgxAgx_dim-1)*(AtgxAgx_dim-1)];
			memset(Q11_right, 0, sizeof(double) * (AtgxAgx_dim-1)*(AtgxAgx_dim-1));
			dger(AtgxAgx_dim-1, AtgxAgx_dim-1, 1/S, iA11A12, 1, A21iA11, 1, Q11_right, AtgxAgx_dim-1);
			
//			iAtgxAgx = zeros(n);
//	        iAtgxAgx(1:n-1,1:n-1) = iA11+ Q11_right;
//	        iAtgxAgx(1:n-1,n) = -iA11A12/S;
//	        iAtgxAgx(n,1:n-1) =  -A21iA11/S;
//	        iAtgxAgx(n,n) = 1/S;
			for(int k = 0 ; k < AtgxAgx_dim-1; k++){
				for(int l = 0 ; l < AtgxAgx_dim-1; l++){
					iAtgxAgx[k*ldAtgxAgx+l] += Q11_right[k*(AtgxAgx_dim-1)+l];
				}
				iAtgxAgx[(AtgxAgx_dim-1)*ldAtgxAgx+k] = -1/S * iA11A12[k];
				iAtgxAgx[k*ldAtgxAgx + AtgxAgx_dim-1] = -1/S * A21iA11[k];
			}
			iAtgxAgx[(AtgxAgx_dim-1)*ldAtgxAgx + (AtgxAgx_dim-1)] = 1/S;
						
//	        xk_1(i_delta) = 0;
			xk_1[i_delta] = 0;
			
			delete [] tempD;
			delete [] tempD2;
			delete [] iA11A12;
			delete [] A21iA11;
			delete [] Q11_right;
			delete [] AtgxAnx;
		}
		
		memset(z_xk, 0, N * sizeof(int));
		for(int k = 0 ; k < gamma_xk_size; k++){
			z_xk[gamma_xk[k]] = -1 * sign(Primal_constrk[gamma_xk[k]]);
		}
		for(int k = 0 ; k < gamma_x_size; k++){
			Primal_constrk[gamma_x[k]] = sign(Primal_constrk[gamma_x[k]]) * epsilon;
		}
		// END OF METHOD START FREEING
		delete [] del_x;
		delete [] Asupported;
	}
	
	nIter = iter;
	
	
	delete [] gamma_x;
	delete [] z_x;
	delete [] xk_temp;
	delete [] del_x_vec;
	delete [] pk_temp;
	delete [] dk;
	delete [] x_k;
	delete [] z_xk;
	delete [] gamma_xk;
	delete [] Primal_constrk;
	delete [] n_zx;
	delete [] AtgxAgx;
	delete [] iAtgxAgx;
	delete [] Agdelx;
	
	x = xk_1;
	return x;
}

void update_primal(int &i_delta, double &delta,int &out_x){
	bool *gamma_lc = new bool[N];
	double *delta_pos = new double[N];
	int delta_pos_len;
	memset(gamma_lc, true, sizeof(bool)*N);
	double min_value;
	int min_index;
	double temp;
	
	int *delta1_pos_ind, *delta2_pos_ind, *delta3_pos_ind;
	
	//	gamma_lc = setdiff(1:N, union(gamma_x, out_x));
	for(int k = 0 ; k < gamma_x_size; k++){
		gamma_lc[gamma_x[k]] = false;
	}
	if (out_x > -1 &&  out_x < N)
		gamma_lc[out_x] = false;
	
//	  delta1_constr = (epsilon-pk_temp(gamma_lc))./(1+dk(gamma_lc));
//    delta1_pos_ind = find(delta1_constr>0);
//    delta1_pos = delta1_constr(delta1_pos_ind);
//    [delta1 i_delta1] = min(delta1_pos);
//    if isempty(delta1)
//        delta1 = inf;
//    end.
	double delta1;
	int i_delta1;
	delta_pos_len = 0;
	min_index = 0;
	delta1_pos_ind = new int[N - gamma_x_size];
	
	for(int k = 0; k < N; k++){
		if (gamma_lc[k] == true){
			temp = (epsilon - pk_temp[k])/(1+dk[k]);			
			if ( temp > 0 ){
				delta1_pos_ind[delta_pos_len] = k;
				delta_pos[delta_pos_len] = temp;
				delta_pos_len++;
			}
		}
	}

	min_value = delta_pos[min_index];
	for(int k = 1 ; k < delta_pos_len; k++){
		if (delta_pos[k] < min_value){
			min_index = k;
			min_value = delta_pos[min_index];
		}
	}
	delta1   = min_value;
	i_delta1 = min_index;
	
//	delta2_constr = (epsilon+pk_temp(gamma_lc))./(1-dk(gamma_lc));
//	delta2_pos_ind = find(delta2_constr>0);
//	delta2_pos = delta2_constr(delta2_pos_ind);
//	[delta2 i_delta2] = min(delta2_pos);
//	if isempty(delta2)
//	    delta2 = inf;
//	end
	double delta2;
	int i_delta2;
	delta_pos_len = 0;
	min_index = 0;
	delta2_pos_ind = new int[N - gamma_x_size];
	
	for(int k = 0; k < N; k++){
		if (gamma_lc[k] == true){
			temp = (epsilon + pk_temp[k])/(1-dk[k]);
			if ( temp > 0 ){
				delta2_pos_ind[delta_pos_len] = k;
				delta_pos[delta_pos_len] = temp;
				delta_pos_len++;
			}
		}
	}
	
	min_value = delta_pos[min_index];
	for(int k = 1 ; k < delta_pos_len; k++){
		if (delta_pos[k] < min_value){
			min_index = k;
			min_value = delta_pos[min_index];
		}
	}
	delta2   = min_value;
	i_delta2 = min_index;
		
//	if delta1>delta2
//	    delta = delta2;
//	    i_delta = gamma_lc(delta2_pos_ind(i_delta2));
//	else
//	    delta = delta1;
//	    i_delta = gamma_lc(delta1_pos_ind(i_delta1));
//	end
	if (delta1 > delta2){
		delta   = delta2;
		i_delta = delta2_pos_ind[i_delta2];
	}else{
		delta   = delta1;
		i_delta = delta1_pos_ind[i_delta1];
	}
	
//	delta3_constr = (-xk_temp(gamma_x)./del_x_vec(gamma_x));
//	delta3_pos_index = find(delta3_constr>0);
//  delta3_pos = delta3_constr(delta3_pos_index);
//	[delta3 i_delta3] = min(delta3_pos);
//	out_x_index = gamma_x(delta3_pos_index(i_delta3));	
	double delta3;
	int i_delta3;
	delta_pos_len = 0;
	min_index = 0;
	delta3_pos_ind = new int[gamma_x_size];
	for(int k = 0 ; k < gamma_x_size; k++){
		temp = -1 * xk_temp[gamma_x[k]] / del_x_vec[gamma_x[k]];
		if (temp > 0){
			delta3_pos_ind[delta_pos_len] = k;
			delta_pos[delta_pos_len] = temp;
			delta_pos_len++;
//			delta3_pos_ind[delta_pos_len] = gamma_x[k];
//			delta_pos[delta_pos_len] = temp;
//			delta_pos_len++;
		}
	}
	
	if (delta_pos_len > 0){
		min_value = delta_pos[min_index];
		for(int k = 1 ; k < delta_pos_len; k++){
			if (delta_pos[k] < min_value){
				min_index = k;
				min_value = delta_pos[min_index];
			}
		}
		delta3   = min_value;
		i_delta3 = min_index;
	}else{
		delta3 = -1;
		i_delta3 = -1;
	}
	
//	printf("delta_pos_len: %d\n", delta_pos_len);
	
//	out_x = [];
//	if ~isempty(delta3) && (delta3 > 0) && (delta3 <= delta)
//	    delta = delta3;
//	    out_x = out_x_index;
//	end
	out_x = -1;
	if (delta3 > 0 && delta3 <= delta){
		delta = delta3;
		out_x = gamma_x[delta3_pos_ind[i_delta3]];
	}
	
//	printf("\t delta1 %f delta2 %f delta3 %f\n", delta1, delta2, delta3);
//	xk_1 = xk_temp+delta*del_x_vec;
//	xk_1(out_x) = 0;
//	daxpy(N, delta, del_x_vec, 1, xk_temp, 1);
//	if (out_x != -1)
//		xk_1[out_x] = 0;

//	wrong_sign = find(sign(xk_1(gamma_x)).*z_x(gamma_x)==-1);
//	int *wrong_sign = new int[gamma_x_size];
//	int wrong_sign_len = 0;
//	for(int k = 0 ; k < gamma_x_size){
//		if ( sign(xk_1[gamma_x]) * z_x[gamma_x[k]] == -1 ){
//			wrong_sign[wrong_sign_len] = k;
//			wrong_sign_len++;
//		}
//	}
		
//	if ~isempty(gamma_x(wrong_sign))
//	    delta = 0;
//	    % can also choose specific element which became non-zero first but all
//	    % that matters here is AtA(gx,gl) doesn't become singular.
//	    [val_wrong_x ind_wrong_x] =  sort(abs(del_x_vec(gamma_x(wrong_sign))),'descend');
//	    out_x = gamma_x(wrong_sign(ind_wrong_x));
//	end
// 
	delete [] gamma_lc;
	delete [] delta_pos;
	delete [] delta1_pos_ind;
	delete [] delta2_pos_ind;
	delete [] delta3_pos_ind;
}
