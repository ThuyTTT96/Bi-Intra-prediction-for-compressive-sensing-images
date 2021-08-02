
%mex -DCOMPILE_MEX -g -largeArrayDims -cxx SolveDALM.cc -lmwblas -lmwlapack
%! mv SolveDALM.`mexext` SolveDALM_mex.`mexext`
%! cp ../regular_CUDA/SolveDALM.`mexext` SolveDALM_cuda.`mexext`
%unix(['mv SolveDALM.' mexext  ' SolveDALM_mex.' mexext]);
%unix(['cp ../regular_CUDA/SolveDALM.' mexext  ' SolveDALM_cuda.' mexext]);

% if ~exist('matlabdir', 'var')
%     matlabdir = '/Applications/MATLAB_R2009b.app';
% end
% 
% cd ../regular_CUDA/
% unix(['nvcc -DCOMPILE_MEX -c -m 64 SolveDALM.cu -Xcompiler -fPIC -I ' matlabdir '/extern/include -I /usr/local/cuda/include/']);
% 
% if exist('/usr/local/cuda/lib64', 'dir')
%     mex -cxx -L/usr/local/cuda/lib -lcuda -lcudart -lcublas -lmwblas -lmwlapack SolveDALM.o;
% else
%     mex -cxx -L/usr/local/cuda/lib -lcuda -lcudart -lcublas -lmwblas -lmwlapack SolveDALM.o;
% end
% cd ../regular
% unix(['mv ../regular_CUDA/SolveDALM.' mexext  ' SolveDALM_cuda.' mexext]);


rand ('seed',0);
randn ('seed',0);

m = 200;
n = 2000;
nnz = 38;
nne = 10;

A = randn (m,n);
x = zeros (n,1);

while sum(x~=0) < nnz,
	k = nnz-sum(x~=0);
	s = 1+floor(n*rand(k,1));
	x(s) = randn(k,1);
end

e = zeros (m,1);
while sum(e~=0) < nne,
	k = nne-sum(e~=0);
	s = 1+floor(m*rand(k,1));
	e(s) = randn(k,1);
end

b = A*x+e;

nu = 10;
tol = 1e-3;
stopCrit = 5;

start = tic;
[xm nIterm] = SolveDALM (A, b, 'lambda', nu, 'tolerance', tol, 'stoppingCriterion', stopCrit, 'groundtruth', x);
['elapsed MATLAB: ' num2str(toc(start))]
start = tic;
[xx nIterx] = SolveDALM_mex (A, b, nu, tol, stopCrit, x);
['elapsed MEX   : ' num2str(toc(start))]

r_mat = norm(A*xm - b);
r_mex = norm(A*xx - b);


close all;


figure,
plot (x, '--bx');
%hold on;
%plot (e, '--g.');
grid on;
title('Truth');

figure,
grid on;
plot (xx, '--bx');
grid on;
title(sprintf('Mex. residual = %f', r_mex));

figure,
plot (xm, '--bx');
title(sprintf('Matlab. residual = %f', r_mat));
grid on;

fprintf('||Xmex - Xmatlab||2 = '); disp(norm(xx-xm));
fprintf('||Xmex - Xmatlab||1 = '); disp(sum(abs(xx-xm)));

fprintf('nIter mex = '); disp(nIterx);
fprintf('nIter matlab = '); disp(nIterm);