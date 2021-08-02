% MEX/C code for PALM test code

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


% this script just runs through different sizes of the A matrix and compares the accuracy and 
% amount of time it takes to solve each problem to the matlab implementation

ratiodn = [10 2];
ks = 0.1;

% elements in tols match up to tol_ints.  
%  it runs tol(1) with tol_ints(1), tol(2) with tol_ints(2)...
tols = [1e-3 1e-4 1e-5 1e-6];
tol_ints = [1e-3 1e-4 1e-5 1e-6];

maxIter = 5000;
maxIter_apg = 50;
iter = 0;
iterations = 10;
dimsPowers = [500:100:900 1000:1000:10000];

stopCrit = 5;

for aa = 1:numel(tols)
    disp(['tolerance ' num2str(tols(aa))]);
    %    for bb = 1:numel(tol_ints)
    bb = aa;
    
    tol = tols(aa);
    tol_int = tol_ints(bb);
    
    for i=1:numel(dimsPowers)
        d = dimsPowers(i);
        n = ceil((d / ratiodn(1)) * ratiodn(2));
        k = ceil(n * ks);
        disp(['dimension of ' num2str(d) ' by ' num2str(n)]);
        
        RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));
        tempM = randn(d, n);
        %matrixNorm = tempM.'*tempM;
        %matrixNorm = sqrt(diag(matrixNorm)).';
        %tempM = tempM./repmat(matrixNorm, [d,1]);
        A = tempM;
        
        iter = 0;
        for j = 1:iterations
            iter = iter + 1;
            %disp(['iteration ' num2str(j)]);
            sparseSupport = randperm(n);
            x0=zeros(n,1);
            x0(sparseSupport(1:k))=randn(1,k);
            
            sparseSupport = randperm(d);
            yk=zeros(d,1);
            yk(sparseSupport(1:k)) = randn(1,k);
            
            x0 = x0 / norm(x0);
            yk = yk / norm(yk);
            
            y0 = A*x0 + yk;
            
            b1 = y0;
            
            %G = A'*A ;
            %opts.disp = 0;
            %tau = eigs(double(G),1,'lm',opts);
            
            a=tic;
            [xm em nIterm nIter_in_total] = SolvePALM_mex(b1, A, tol, tol_int, maxIter, maxIter_apg);
            [xm1 em1 nIterm1 nIter_in_total1] = norm_x_e_primal(b1, A, tol, tol_int, maxIter, maxIter_apg);
            
            disp(['norm(xm-xm1): ' num2str(norm(xm-xm1))]);
            disp(['norm(em-em1): ' num2str(norm(em-em1))]);
            disp(['nIterm-nIterm1: ' num2str(nIterm-nIterm1)]);
            disp(['nIterT-nIterT1: ' num2str(nIter_in_total-nIter_in_total1)]);
            
        end
        
    end
end