/*% C code for L1-min test code

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

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "portable_blas_wrapper.h"
#include "timing.h"
#include "SolveHomotopy.h"

using namespace std;

int main(int argc, char **argv) {
        TIMINGINIT
			int alg = 1;
	int m = 1024;
	int n = 32*10;
	double sparsity = 0.1;
	double tol = 0.01;
	int iter, iter_total;
	
	int seed = 100;
	
	srand(seed);
	
	double lambda = 0.01;
	int maxIter = 5000;
	double *b  = new double[m];
	double *yk = new double[m];
	double *A  = new double[m*n];
	double *x  = new double[n];
	double *xG = new double[n];
	double *e  = new double[m];
	double diffX;
	double diffE;
	double normX;
	double normXG;
	
	double normB, normX0, normA;
	normB  = 0;
	normX0 = 0;
	
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

	dgemv('N', m, n, 1, A, m, xG, 1, 0, b, 1);
	
	for(int i = 0; i < m; i++){
		normB += b[i]*b[i];
	}
		
	cout << "Allocating a solver instance..." << endl;
	
    cout << "ALGORITHM HOMOTOPY min ||x||1 " << endl;

	cout << "Staring minimization... " << endl;
	TIC
	SolveHomotopy(x, iter, b, A, tol, lambda, maxIter, m, n);
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
	cout << "Instance " << 1 << " had " << iter << " iteractions and had norm1(xG-x,1) of " << diffX << " and norm1(x) of " << normX<< " and norm1(xG) of " << normXG << endl;
	
	cout << "Code took " << timeDifference << " sec per instance to run!" << endl;
	cout << "Freeing Memory...\n";
	delete [] A;
	delete [] b;
	delete [] x;
	delete [] xG;
	delete [] e;
	cout << "...Done!\n";
}
