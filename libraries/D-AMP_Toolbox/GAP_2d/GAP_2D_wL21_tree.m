function [theta model]  = GAP_2D_wL21_tree(y0, P, Phi, row,col, block, domain, m_star,stopc,acc)
% Xin Yuan, try to repeat Xuejun's GAP 2D
% initial date: 04/30/2013

% Input: y, the 1D measurement
% Input: P, the P = Phi'*inv(Phi*Phi'), saved matrix for update
% m_star, the index number used in GAP for the non_zero and zero coeffients
% decomposition

% Note: here we use the block size by the L_2,1 norm

[xdim ydim] = size(P);
row_block = block.row;
col_block = block.col;


wave_level = 3;
nlogx = log2(row);
nlogy = log2(col);

Newblock.x = 2.^(wave_level:1:(nlogx-1));
Newblock.y = 2.^(wave_level:1:(nlogy-1));

Newblock.x = [2.^(wave_level) Newblock.x];
Newblock.y = [2.^(wave_level) Newblock.y];

if nargin<10
    acc=0;
end

if nargin<9
    stopc.iternum = 10^2;
    stopc.err = 10^-5;
end

if nargin<8
    m_star = ceil(ydim/(row_block*col_block));
end

if nargin<7
    domain = 'dct';
end

iternum = stopc.iternum;
iterr = stopc.err;
verbose = 'on';



theta = zeros(xdim,1); % time domain signal
w = zeros(xdim,1);   % Inital value, frequency domain sparse
y = zeros(ydim,1);


for iter_i = 1:(iternum)
    theta_old = theta; % Used to check the stop criterion
    if(acc)
        y = y + (y0 - Phi*theta);
    else
        y = y0;
    end
    theta_temp = theta+ P*(y-Phi*theta);
    
    if strcmp(domain,'dct')
        w = myfastdct2(theta_temp,row,col);
    else
        w = myfastdwt2_level(theta_temp,row,col,wave_level);
    end

    % unlike the L_1 norm, here we need to sort the group L_2,1 norm
    % Now we compute the L2,1 norm for each group
    w_temp = reshape(w,[row col]);
    % The followign things is to get weight
    w_temp_vec = im2col(w_temp, [row_block col_block],'distinct');
    
    %weight = my_get_weight(row,col,block);
    weight = my_get_weight_wavelettree(row, col,Newblock);
    w_temp_weight = w_temp.*weight;
    
    w_temp_weight_vec = im2col(w_temp_weight, [row_block col_block],'distinct');
    w_temp_vec_l21 = sqrt(sum(w_temp_weight_vec.^2)); 
    
    % Now we sort the w_temp_vec
    w_sort = sort(w_temp_vec_l21,'descend');
    lambda = w_sort(m_star+1);
     w_temp_vec_l21(find(w_temp_vec_l21<=1e-8)) = lambda;
    w_temp_vec_l21_minus = (max(w_temp_vec_l21-lambda,0))./w_temp_vec_l21;
    w_temp_vec_post = w_temp_vec.*repmat(w_temp_vec_l21_minus,[row_block*col_block ,1]);
    w_temp_post = col2im(w_temp_vec_post,[row_block col_block],[row col],'distinct');
    %w_temp = sign(w).*max(abs(w)-lambda,0);

    if strcmp(domain,'dct')
        theta = myfastidct2(w_temp_post(:),row,col);
    else
        theta = myfastidwt2_level(w_temp_post(:),row,col,wave_level);
    end
    
    g(iter_i) = sqrt(sum((theta_old-theta).^2));

    if(mod(iter_i,20)==1)
    if strcmp(verbose,'on')
        fprintf('GAP_weighted_L2,1 %d of %d: distance = %1.10f,\n',iter_i,iternum,g(iter_i));
    end
    end
    if iter_i>2 & g(end-2:end)<(iterr)      
        break;
    end
    
end
model.w = w;
model.theta = theta;
model.g = g;


end