clear 
close all
clc

%format short

% GAP for mosaic and CACTI reconstruction
% Xin Yuan; ECE, Duke
% xin.yuan@duke.edu
% initial date: 2013-06-04

fname = 'Traffic';
load(fname);
Xtst = X(1:256,1:256,:);  % original video

CodeFrame = 8;  % How many coded images, measurement
spec_channel = 4;  % How many frames are collasped into one measurement
[Row, Col, Frame] = size(Xtst);

ColT = spec_channel;

% generate binary code
 Phi = binornd(1,0.5,[Row, Col,ColT]);
 Phi(:,:,1) = binornd(1,0.5,[Row, Col]);
 shift = 1;  % we use shift code
 for t=2:ColT
     % Phi(:,1+(t-1)*shift:Col,t) = Phi(:,1+(t-2)*shift:Col-shift,t-1);  % mask is moving from  left to right
     Phi(1+(t-1)*shift:Row,:,t) = Phi(1+(t-2)*shift:Row-shift,:,t-1);    % mask is moving from  up to down
 end

 % generate measurement
 y = zeros(Row,Col,CodeFrame);
for t = 1:CodeFrame
   y(:,:,t) = sum(Xtst(:,:,(t-1)*ColT+(1:ColT)).*Phi,3);
end

stopc.iternum = 30;
stopc.err = 10^-5;
acc = 2;
ydim = Row*Col;


spbasis.space = 'wavelet';  % transform for space, 'wavelet' or 'dct'
spbasis.time  = 'dct';  % transform for spectrum, 'wavelet' or 'dct', dct is always used. If we use wavelet, T need to be the power of 2. we can use no, means no transformation
%spbasis.spectrum  = 'no';  % Here we use no, means no transfromation in spectrum

weighttype.space = 'tree';   % Here we can select:  'tree' or 'block'
weighttype.time = 'block';   % Here we can select:  'tree' or 'block', if we use tree, T need to be the power of 2


weight_base.type = 'cos'; % here we can select 'exp' or 'cos'
if strcmp(weight_base.type,'exp')
    weight_base.space = 2;   % This weight is the base of exponential decay. should be large than 1 [1 2] is always used
    weight_base.time = 2;
end



% The block size for group
block.row = 2;
block.col = 2;
block.T = 2;

m_star = ceil(ydim/(block.row*block.col*block.T));

    for k=1:CodeFrame        
        fprintf('----- Reconstructing frame-block %d of %d\n',k,CodeFrame);
     
        y_use = y(:,:,k);
        
        theta_wL21 = GAP_3D_wL21_grayscale(y_use,Phi, Row,Col,spec_channel, block,spbasis, m_star,stopc,acc,weight_base,weighttype); 
  
        img(:,:,(k-1)*spec_channel+(1:spec_channel)) = theta_wL21;
      
    end

    % show result
    figure;
 for t=1:spec_channel
      PSNR(t) = SS_PSNR_3D(Xtst(:, :,t),img(:,:,t));
     imshow([Xtst(:, :,t) img(:,:,t)]); title(num2str(t)); pause(0.5); 
 end
 
     % plot PSNR
 figure; plot(PSNR);

