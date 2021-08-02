/*% MEX/C code for DALM l1-minimization

% Copyright ©2010. The Regents of the University of California (Regents).
% All Rights Reserved. Contact The Office of Technology Licensing,
% UC Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA 94720-1620,
% (510) 643-7201, for commercial licensing opportunities.

% Created by Victor Shia, Allen Y. Yang, Department of EECS, University of California,
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

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "portable_blas_wrapper.h"
#include "timing.h"
#include "SolveDALM_fast.h"

#define double double

using namespace std;

int main(int argc, char **argv) {
        TIMINGINIT
			int alg = 1;
	int m = 1024;
	int n = 32*10;
	int maxIter = 5000;
	double sparsity = 0.1;
	double tol = 0.001;
	int iter, iter_total;
	
	int seed = 100;
	
	srand(seed);
	
	double lambda = 0.01;
	double *b = new double[m];
	double *yk = new double[m];
	double *A = new double[m*n];
	double *x = new double[n];
	double *xG = new double[n];
	double *e;
	double diffX;
	double diffE;
	double normX;
	double normXG;
	
	double normB, normX0, normA, normYk;
	normB  = 0;
	normX0 = 0;
	normYk = 0;	
	
	for(int i = 0 ; i < m; i++){
		normA  = 0;
		for(int j = 0 ; j < n; j++){
			A[j*m + i] = (double)1.0 * (rand() - rand())/100 ;
			normA += A[j*m + i]*A[j*m + i];
		}
		normA = sqrt(normA);
		for(int j = 0 ; j < n; j++){
			A[j*m + i] = (double)(A[j*m + i] / normA);
		}
	}
	normA  = 0;
	for(int i = 0 ; i < m; i++){
		for(int j = 0 ; j < n; j++){
			normA += A[j*m + i]*A[j*m + i];
		}
	}
			
	normX0 = 0;
	for(int j = 0 ; j < n; j++){
		xG[j] = 0;
	}
	for(int j = 0 ; j < (int)(n * sparsity); j++){
		xG[(int)(rand()*n) % n] = (double)1.0 * (rand() - rand()) ;
	}
	for(int j = 0; j < n; j++){
		normX0 += xG[j]*xG[j];
	}
	normX0 = sqrt(normX0);
	for(int j = 0 ; j < n; j++){
		xG[j] = (double)(xG[j]/normX0);
	}
	normX0 = 0;
	for(int j = 0; j < n; j++){
		normX0 += xG[j]*xG[j];
	}
	
	normYk = 0;
	for(int j = 0 ; j < m; j++){
		yk[j] = 0;
	}
	for(int j = 0 ; j < (int)(m * sparsity); j++){
		yk[(int)(rand()*m) % m] = (double)1.0 * (rand() - rand()) ;
	}
	for(int j = 0; j < m; j++){
		normYk += yk[j]*yk[j];
	}
	
	normYk = sqrt(normYk);
	
	for(int j = 0 ; j < m; j++){
		//e[j] = 0;
		yk[j] = (double)(yk[j]/normYk);
	}
	
	normYk = 0;
	for(int j = 0; j < m; j++){
		normYk += yk[j]*yk[j];
	}
	
	dgemv('N', m, n, 1, A, m, xG, 1, 0, b, 1);

	if (alg == 3){
	for(int i = 0 ; i < m ; i++){
		b[i] += yk[i];
	}
	}
	
	for(int i = 0; i < m; i++){
		normB += b[i]*b[i];
	}
		
	cout << "Allocating a solver instance..." << endl;
	
    cout << "ALGORITHM DALM_fast min ||x||1" << endl;

	cout << "Staring minimization... " << endl;
	TIC
	SolveDALM_fast(x, iter, b, A, lambda, tol, maxIter, m, n, 5, xG);
	TOC
	
	// check accuracy
	diffX = 0;
	diffE = 0;
	normX = 0;
	normXG = 0;	
	for(int j = 0 ; j < n; j++){
		normX += fabs(x[j]);
		normXG += fabs(xG[j]);
		diffX += fabs(xG[j] - x[j]);
	}
	if (alg == 3){
	  for(int j = 0 ; j < m; j++){
	    diffE += fabs(e[j] - yk[j]);
	  }
	}
	cout << "Instance " << 1 << " had " << iter << " iteractions and had norm1(xG-x,1) of " << diffX << " and norm1(x) of " << normX<< " and norm1(xG) of " << normXG << " and norm1(yk-e) of " << diffE << "\n";
	
	cout << "Code took " << timeDifference << " sec per instance to run!" << endl;
	cout << "Freeing Memory...\n";
	delete [] A;
	delete [] b;
	delete [] yk;
	delete [] x;
	delete [] xG;
	cout << "...Done!\n";
}
