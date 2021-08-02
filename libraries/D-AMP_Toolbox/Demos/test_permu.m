% this code is used for test the permutation
filename='barbara.png';
imsize=256;

ImIn=double(imread(filename));
x_0=imresize(ImIn,imsize/size(ImIn,1));
[height, width]=size(x_0);
load('perm64k.mat');

[B,Q] = sort(P64k);
y_all = fwht(x_0(P64k));
x_0_re = fwht(y_all);
x_0_re = x_0_re(Q);


%figure; imagesc([x_0, reshape(x_0_re, [height, width])]);

figure; imagesc(reshape(x_0_re, [height, width]))