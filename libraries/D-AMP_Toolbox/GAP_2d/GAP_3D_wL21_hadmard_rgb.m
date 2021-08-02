function theta_3d = GAP_3D_wL21_hadmard_rgb(y,Row,Col,T, block,spbasis, ydim,stopc,acc,weight_base,weighttype)
    % Xin Yuan, duke ECE
    % xin.yuan@duke.edu
    % initial date: 04/03/2015
   % Input parameters:
   % ---- y, the measuerment, compressed measurements for each channel
   % -----Phi, the code obtained from the moving mask,
   % ---- Row, Col, T are the videosize to be reconstructed
   % ---- Block, the 3D block size, a structure with row col T block
   % ---- spbasis, transform basis, dct or wavelet, space and time can have
   % different basis
   % ---- m_star, the m_star in the GAP paper
   % ---- stopc, the stop criterion to stop the itration, iternum or error
   % ---- weight_base, if we use the exponential decay weight, we need to
   % select the base, time and space can have different base
   % ---- weighttype.space   can be tree or block
   % --  weighttype.time = 'block';   % Here we can select:  'tree' or 'block', if we use tree, T need to be the power of 2

    
row = Row;
col = Col;


row_block = block.row;
col_block = block.col;
T_block = block.T;
%ydim = m_star;
m_star = ceil(3*ydim/(row_block*col_block*T_block));

if nargin<11
    weighttype.space = 'tree';  % will affect the weight selection
    weighttype.time = 'block';
end

if nargin<10
    weight_base='cos';  % will affect the weight selection
end

if nargin<9
    acc=0;
end

if nargin<8
    stopc.iternum = 10^2;
    stopc.err = 10^-5;
end

if nargin<7
    m_star = ceil(ydim/(row_block*col_block*T_block));
end

if nargin<6
    spbasis.space = 'wavelet';
    spbasis.time = 'dct';
end


iternum = stopc.iternum;
iterr = stopc.err;
verbose = 'on';

%Phi_sum = sum(Phi.^2,3);



if strcmp(weighttype.time,'block')
    block_num_T = ceil(T/T_block);
    T_new = T_block*block_num_T;
    if strcmp(weighttype.space,'block')
        if strcmp(weight_base.type,'cos')
            weight = my_get_weight3d(row,col,T_new,block);
        else
            weight = my_get_weight3d_exp(row,col,T_new,block,weight_base);
        end
    else
            wave_level_space = 3;
            nlogx = log2(row);
            nlogy = log2(col);
            Newblock.row = 2.^(wave_level_space:1:(nlogx-1));
            Newblock.col = 2.^(wave_level_space:1:(nlogy-1));

            Newblock.row = [2.^(wave_level_space) Newblock.row];
            Newblock.col = [2.^(wave_level_space) Newblock.col];
            Newblock.T = T_block;
        if strcmp(weight_base.type,'cos')
            weight =  my_get_weight3d_tree(row,col,T_new,Newblock);
        else
            weight = my_get_weight3d_tree_exp(row,col,T_new,Newblock,weight_base);
        end
    end
else
    wave_level_space = 3;
    wave_level_time = 1;
    nlogx = log2(row);
    nlogy = log2(col);
    nlogt = ceil(log2(T));
    T_new = 2^nlogt;
    block_num_T = ceil(T_new/T_block);
    
    Newblock.x = 2.^(wave_level_space:1:(nlogx-1));
    Newblock.y = 2.^(wave_level_space:1:(nlogy-1));
    Newblock.T = 2.^(wave_level_time:1:(nlogt-1));

    Newblock.x = [2.^(wave_level_space) Newblock.x];
    Newblock.y = [2.^(wave_level_space) Newblock.y];
    Newblock.T = [2.^(wave_level_time) Newblock.T];

    weight =  my_get_weight3d_3dtree_exp(row,col,T_new,Newblock,weight_base); 

end

index = reshape(1:row*col*T_new,row,col,T_new);

index_2d1 = im2col(index(:,:,1),[row_block col_block],'distinct');
index_block2 = repmat(index_2d1, [T_block ,1]) + kron([0: (T_block-1)]', row*col*ones(size(index_2d1)));

grid_idx= repmat(index_block2 , [1, block_num_T]) + kron([0: (block_num_T-1)], row*col*T_block*ones(size(index_block2)));
clear index_2d1  index_block2


% get transform matrix
 if strcmp(spbasis.space,'dct')
  T_row = dct(eye(row));
  T_col = dct(eye(col));  
else
 level = 3;
 qmf   = MakeONFilter('Daubechies',8); 
 sig_level_row = log2(row); 
 sig_level_col = log2(col); 
 T_row = get_waveletMatrix(qmf,sig_level_row,level,level);
 T_col = get_waveletMatrix(qmf,sig_level_col,level,level);
end

if strcmp(spbasis.time,'dct')
    T_t = dct(eye(T));
else
    T = 2^(ceil(log2(T)));

   level = 1;
   qmf   = MakeONFilter('Haar',4);  % Here we use the 'Haar' wavelet and the second parameter now no meaning
   sig_level_t = ceil(log2(T));
   T_t = get_waveletMatrix(qmf,sig_level_t,level,level);

end


theta_3d = zeros(row*col,T);




for iter_i = 1:(iternum)
    theta_3d_old = theta_3d; % Used to check the stop criterion
      tic;
     %theta_3d = reshape(theta,[row, col,T]);
       y0_2d  = y;
         
        if(acc)
            for nr = 1:T
            ytemp = y0_2d(:,nr) - getPhitheta(theta_3d(:,nr),ydim);
            y_2d(:,nr) = y0_2d(:,nr) + ytemp;
            end
        else
            y_2d = y0_2d;
        end
        
        for nr = 1:T
        theta_temp(:,nr) = theta_3d(:,nr) + getPinvPhitheta(y_2d(:,nr),theta_3d(:,nr),ydim);
        end
        theta_temp_3d = reshape(theta_temp,[row, col, T]);
        w_3d = myfasttransform3_givenT(theta_temp_3d,T_row,T_col,T_t);  % Note here theta_temp and w are all 3D

       % weight = my_get_weight3d(row,col,T,block);
         if(T_new~=T)
            w_3d(:,:,(T+1):T_new) = repmat(squeeze(w_3d(:,:,T)),[1,1,(T_new-T)]);
         end
         
         w_3d_weight = w_3d.*weight;
         block_w_3d = w_3d_weight(grid_idx);
         block_sum = sqrt(sum(block_w_3d.^2,1));
         
         sum_sort = sort(block_sum,'descend');
        lambda = sum_sort(m_star+1);  
        block_sum(find(block_sum<=1e-8)) = lambda;
        block_sum_minus = max(block_sum-lambda,0)./block_sum;

        block_sum_minus_vec = repmat(block_sum_minus,[size(block_w_3d,1) ,1]);

        w_3d = w_3d(grid_idx);
        w_3d_post = w_3d.*block_sum_minus_vec;
        w_3d_post1(grid_idx) = w_3d_post;
        % Now we need to reshape back
        w_temp_post = reshape(w_3d_post1,[row, col, T_new]);
        if(T_new>T)
            w_temp_post = w_temp_post(:,:,1:T);
        end


        theta_3d = reshape(myfastinvtransform3_givenT(w_temp_post,T_row,T_col,T_t),[row*col, T_new]);

      

        for nr = 1:T
        L2_err(nr) = norm(y0_2d(:,nr) - getPhitheta(theta_3d(:,nr),ydim));
        end
      
   
         ittime = toc;
          if(iter_i==1)
              disp(['One iteration time: ' num2str(ittime)]);
          end
        meanL2_err = norm(L2_err);
        clear L2_err;
        g(iter_i) = sqrt(sum(sum((theta_3d_old-theta_3d).^2)));
        
        if(mod(iter_i,20)==1)
        if strcmp(verbose,'on')
            fprintf('GAP3d_L2,1 %d of %d: distance = %1.10f, L2 error = %1.10f\n',iter_i,iternum,g(iter_i),meanL2_err);
        end
        end
        if iter_i>2 & g(end-2:end)<(iterr)      
            break;
        end
    
end


end