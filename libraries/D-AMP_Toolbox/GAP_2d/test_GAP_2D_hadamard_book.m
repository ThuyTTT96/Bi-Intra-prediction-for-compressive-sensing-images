close all
clear 
clc


% GAP for 2D signal
% Xin Yuan 
% initial Date: 04/30/2013
% Model: y=Phi x, where x is sparse

filename = 'Books';
X2D = im2double(rgb2gray(imread([filename '.jpg'])));
figure; imshow(X2D);
Row = 128;
Col  =128;

X2D = imresize(X2D,[Row, Col]);


[row col] = size(X2D);

x = X2D(:);
xdim = size(x,1);

%y_cand = 20:8:64;

%y_cand = 30;

%for nc = 1:length(y_cand)
CSratio = 0.25;
ydim = round(CSratio*xdim);
% Phi = randn(ydim,xdim);
% Phi = Phi./repmat(sqrt(sum(Phi.^2,1)),[ydim,1]);	
%Phi_all = fwht(eye(xdim));
%Phi = Phi_all(1:ydim,:);
load('Hardmard16k.mat');%load('Perm16K');
%Phi_all = Phi_all(P,:);
Phi = Phi_all(1:ydim,:);
% now we do permutation

%Phi = Phi(:,P);
clear Phi_all


%Phi = randn(ydim, xdim);
y = Phi*x;  % generate the measurement
%y_all = fwht(x);
%y_check = y_all(1:ydim,:);


% The following we need to call the GAP
% First, we save the P
%P = Phi'*inv(Phi*Phi');
Pinv = Phi';
%%
stopc.iternum = 80;
stopc.err = 10^-5;
acc = 1;

% The L1 
m_star_L1 = ydim;

% [theta_L1_dct modelL1_dct] = GAP_2D_L1(y, P, Phi, row,col, 'dct', m_star_L1,stopc,acc);
% [theta_L1_dwt modelL1_dwt] = GAP_2D_L1(y, P, Phi, row,col, 'wavelet', m_star_L1,stopc,acc);

block.row = 2;
block.col = 2;
m_star = ceil(ydim/(block.row*block.col));

% The L_{2,1}
% [theta_L21_dct model_L21_dct] = GAP_2D_L21(y, P, Phi, row,col, block,'dct', m_star,stopc,acc);
% [theta_L21_dwt model_L21_dwt] = GAP_2D_L21(y, P, Phi, row,col, block,'wavelet', m_star,stopc,acc);
% 
% % The weighted L_{2,1}
% [theta_wL21_dct model_wL21_dct] = GAP_2D_wL21(y, P, Phi, row,col, block,'dct', m_star,stopc,acc);
% [theta_wL21_dwt model_wL21_dwt] = GAP_2D_wL21(y, P, Phi, row,col, block,'wavelet', m_star,stopc,acc);
[theta_wL21_dwt_tree model_wL21_dwt_tree] = GAP_2D_wL21_tree(y, Pinv, Phi, row,col, block,'wavelet', m_star,stopc,acc);



% 
% 
% theta_2d_wL21_dct = reshape(theta_wL21_dct,[row, col]);
% PSNR_wL21_dct(nc) = SS_PSNR_3D(X2D,theta_2d_wL21_dct);
% 
% theta_2d_wL21_dwt = reshape(theta_wL21_dwt,[row, col]);
% PSNR_wL21_dwt(nc) = SS_PSNR_3D(X2D,theta_2d_wL21_dwt);

theta_2d_wL21_dwt_tree = reshape(theta_wL21_dwt_tree,[row, col]);
PSNR_wL21_dwt_tree = psnr(X2D,theta_2d_wL21_dwt_tree/max(theta_2d_wL21_dwt_tree(:)));

%figure; imshow([[theta_2d_L1_dct theta_2d_L21_dct theta_2d_wL21_dct]; [theta_2d_L1_dwt theta_2d_L21_dwt theta_2d_wL21_dwt]]); title(['weighted L_{2,1}, L_{2,1} L_1 GAP, CS ratio: ' num2str(CSratio(nc)) ' PSNR: ' num2str(PSNR_wL21_dwt(nc))])
%figure; imshow([w_2d   theta_2d_L1]); title(['L_{2,1} and L_1 GAP, CS ratio: ' num2str(CSratio(nc)) ' PSNR: ' num2str(PSNR(nc))])
%figure; imshow([theta_2d]); title(['weighted L_{2,1} GAP, CS ratio: ' num2str(CSratio(nc)) ' PSNR: ' num2str(PSNR(nc))])


%end

figure; imagesc([theta_2d_wL21_dwt_tree]); colormap gray; title(['weighted L_{2,1},of GAP, CS ratio: ' num2str(CSratio) ' PSNR: ' num2str(PSNR_wL21_dwt_tree)])


% figure;
% 
% plot(CSratio,PSNR_wL21_dwt_tree,'r*-');
% hold on;
% plot(CSratio,PSNR_wL21_dwt,'g>-');
% hold on;
% plot(CSratio,PSNR_wL21_dct,'b-');
% xlabel('CS ratio');
% ylabel('PSNR (dB)');
% title('weighted L_{2,1},GAP, with DCT and wavelet');
% legend( 'wavelet with wavelet-tree structure weights', 'wavelet with block structure weights' ,'DCT with block structure weights');


