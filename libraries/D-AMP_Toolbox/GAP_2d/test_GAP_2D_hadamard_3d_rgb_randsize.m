close all
clear all
clc


% GAP for 2D signal
% Xin Yuan 
% initial Date: 04/30/2013
% Model: y=Phi x, where x is sparse

filename = 'Books';
X2D = im2double((imread([filename '.jpg'])));
%figure; imshow(X2D);
Row = 256;
Col = 200;

X2D = imresize(X2D,[Row, Col]);
figure; imshow(X2D);


[row, col, rgb] = size(X2D);
xdim = row*col;
CSratio = 0.125;
ydim = round(CSratio*xdim);


for nr = 1:rgb
    temp_2d = X2D(:,:,nr);
x = temp_2d(:);
x0 = zeros(2^(ceil(log2(Row*Col))),1);
x0(1:Row*Col) = x;
xdim = size(x,1);



% generate the measurement, use the fwht transfomation for this
% measurement, thus we don't have size limit.
y_all = fwht(x0);
y(:,nr) = y_all(1:ydim,:);
clear y_all
end

% The following we need to call the GAP
% First, we save the P
%P = Phi'*inv(Phi*Phi'); % this is the equation to get Pinv, however, in
%our case, it is equal to Phi'
%Pinv = Phi';
%%
stopc.iternum = 10;
stopc.err = 10^-5;
acc = 1;



spbasis.space = 'wavelet';  % transform for space, 'wavelet' or 'dct'
spbasis.time  = 'dct';  % transform for spectrum, 'wavelet' or 'dct', dct is always used. If we use wavelet, T need to be the power of 2. we can use no, means no transformation
%spbasis.spectrum  = 'no';  % Here we use no, means no transfromation in spectrum

weighttype.space = 'tree';   % Here we can select:  'tree' or 'block'
weighttype.time = 'block';   % Here we can select:  'tree' or 'block', if we use tree, T need to be the power of 2


weight_base.type = 'exp'; % here we can select 'exp' or 'cos'
if strcmp(weight_base.type,'exp')
    weight_base.space = 2;   % This weight is the base of exponential decay. should be large than 1 [1 2] is always used
    weight_base.time = 1;
end



% The block size for group
block.row = 2;
block.col = 2;
block.T = 1;
%theta_3d = GAP_3D_wL21_hadmard_rgb(y,Row,Col,rgb, block,spbasis, ydim,stopc,acc,weight_base,weighttype);
theta_3d = GAP_3D_wL21_hadmard_rgb_randsize(y,Row,Col,rgb, block,spbasis, ydim,stopc,acc,weight_base,weighttype);

for nr = 1:rgb
%temp = reshape(theta_3d(1:Row*Col,nr),[row, col]);
temp = theta_3d(1:Row, 1:Col,nr);
theta_2d_wL21_dwt_tree(:,:,nr) = temp/max(temp(:));
end
%theta_2d_wL21_dwt_tree = theta_3d(1:Row, 1:Col,:);
PSNR_wL21_dwt_tree = psnr(X2D,theta_2d_wL21_dwt_tree/max(theta_2d_wL21_dwt_tree(:)));


%end

figure; imshow([theta_2d_wL21_dwt_tree]); 
title(['weighted L_{2,1},of GAP, CS ratio: ' num2str(CSratio) ' PSNR: ' num2str(PSNR_wL21_dwt_tree)])


