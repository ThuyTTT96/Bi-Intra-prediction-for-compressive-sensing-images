%Demonstrates compressively sampling and D-AMP recovery of an image.
close all
clear all
clc

addpath(genpath('..'));

%Parameters
denoiser1='BM3D';%Available options are NLM, Gauss, Bilateral, BLS-GSM, BM3D, fast-BM3D, and BM3D-SAPCA 
denoiser2='fast-BM3D';
filename='barbara.png';
SamplingRate=.1;
iters=30;
imsize=256;

ImIn=double(imread(filename));
x_0=imresize(ImIn,imsize/size(ImIn,1));
[height, width]=size(x_0);
n=length(x_0(:));
m=round(n*SamplingRate);

y_all = fwht(x_0(:));
y = y_all(1:m,:);
clear y_all

Phi = @(z) A_mywhth_damp_hadamard(z, m);
Phit = @(z) At_mywhth_damp_hadamard(y,n);

%Recover Signal using D-AMP algorithms
x_hat1 = DAMP(y, iters,height, width,denoiser1, Phi, Phit);
x_hat2 = DAMP(y, iters,height, width,denoiser2, Phi, Phit);
% x_hat1 = DAMP(y,iters,height,width,denoiser1,M);
% x_hat2 = DAMP(y,iters,height,width,denoiser2,M);

%D-AMP Recovery Performance
performance1=PSNR(x_0,x_hat1);
performance2=PSNR(x_0,x_hat2);
[num2str(SamplingRate*100),'% Sampling ', denoiser1, '-AMP Reconstruction PSNR=',num2str(performance1)]
[num2str(SamplingRate*100),'% Sampling ', denoiser2, '-AMP Reconstruction PSNR=',num2str(performance2)]

%Plot Recovered Signals
subplot(1,3,1);
imshow(uint8(x_0));title('Original Image');
subplot(1,3,2);
imshow(uint8(x_hat1));title([denoiser1, '-AMP']);
subplot(1,3,3);
imshow(uint8(x_hat2));title([denoiser2, '-AMP']);
