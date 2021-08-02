%Demonstrates compressively sampling and D-AMP recovery of an image.
close all
clear all
clc
addpath(genpath('..'));

%Parameters
denoiser1='BM3D';%Available options are NLM, Gauss, Bilateral, BLS-GSM, BM3D, fast-BM3D, and BM3D-SAPCA 
denoiser2='fast-BM3D';
filename='barbara.png';
SamplingRate=.2;
iters=30;
imsize=256;

ImIn=im2double(imread(filename));
x_0=imresize(ImIn,imsize/size(ImIn,1));
[height, width]=size(x_0);
n=length(x_0(:));
m=round(n*SamplingRate);

row = height;
col = width;
Mea_num = m;
Phi_index = Hadamard_index_zigzag(row*col,Mea_num);
A = @(z) A_wht_ord(z, Phi_index);
At = @(z) At_wht_ord(z, Phi_index, row, col);

% load('perm64k.mat');
% [B,Q] = sort(P64k);
% 
% temp = x_0(:);
% y_all = fwht(temp(Q));
% y = y_all(1:m,1);
y = A(x_0(:));
%clear y_all
%%
% Phi = @(z) A_mywhth_damp_hadamard_perm(z, m,Q);
% Phit = @(z) At_mywhth_damp_hadamard_perm(y,n,P64k);

%Recover Signal using D-AMP algorithms
x_hat1 = DAMP(y, iters,height, width,denoiser1, A, At);
%x_hat2 = DAMP(y, iters,height, width,denoiser2, A, At);
% x_hat1 = DAMP(y,iters,height,width,denoiser1,M);
% x_hat2 = DAMP(y,iters,height,width,denoiser2,M);

%D-AMP Recovery Performance
performance1=psnr(x_0,x_hat1);
%performance2=PSNR(x_0,x_hat2);
[num2str(SamplingRate*100),'% Sampling ', denoiser1, '-AMP Reconstruction PSNR=',num2str(performance1)]
%[num2str(SamplingRate*100),'% Sampling ', denoiser2, '-AMP Reconstruction PSNR=',num2str(performance2)]

%Plot Recovered Signals
subplot(1,3,1);
imshow((x_0));title('Original Image');
subplot(1,3,2);
imshow((x_hat1));title([denoiser1, '-AMP']);
%subplot(1,3,3);
%imshow(uint8(x_hat2));title([denoiser2, '-AMP']);
