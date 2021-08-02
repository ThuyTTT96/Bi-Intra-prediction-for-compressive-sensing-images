function weight = my_get_weight_tree(row,col,block)
% Xin Yuan, get the weight of coefficient for GAP 2D
% initial date: 04/30/2013

% Input: row col, the image size or the coeffient size
% Input: block, the blocksize

% output: weight, will have the size [row col], in order to weight the
% coefficient one by one
% Here we use the cosine weight to weight the coefficient the low frequency
% will have the larger weight

row_block = block.x;
col_block = block.y;

block_num_row = length(row_block);
block_num_col = length(col_block);


weight_row = 1e-5+cos((0:(block_num_row-1))/block_num_row*pi/2).^2;
weight_col = 1e-5+cos((0:(block_num_col-1))/block_num_col*pi/2).^2;

%rowsum = 0;

rowsum = 0; 
for nR =1:block_num_row
    
    colsum = 0;
    for nC = 1:block_num_col
        temp_weight = weight_row(nR) + weight_col(nC);
        weight(rowsum+(1:row_block(nR)), colsum+(1:col_block(nC))) = kron(temp_weight,ones(row_block(nR),col_block(nC)));
        colsum = colsum + col_block(nC);
    end
    rowsum = rowsum + row_block(nR);
end

% weight_row_all = repmat(weight_row', [1, block_num_col]);
% weight_col_all = repmat(weight_col, [block_num_row, 1]);

%weight_block = weight_row_all+weight_col_all;

%weight = kron(weight_block, ones(row_block,col_block));
end