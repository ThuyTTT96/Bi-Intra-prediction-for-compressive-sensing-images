function weight = my_get_weight(row,col,block)
% Xin Yuan, get the weight of coefficient for GAP 2D
% initial date: 04/30/2013

% Input: row col, the image size or the coeffient size
% Input: block, the blocksize

% output: weight, will have the size [row col], in order to weight the
% coefficient one by one
% Here we use the cosine weight to weight the coefficient the low frequency
% will have the larger weight

row_block = block.row;
col_block = block.col;

block_num_row = row/row_block;
block_num_col = col/col_block;


weight_row = 1e-5+cos((0:(block_num_row-1))/block_num_row*pi/2).^2;
weight_col = 1e-5+cos((0:(block_num_col-1))/block_num_col*pi/2).^2;

weight_row_all = repmat(weight_row', [1, block_num_col]);
weight_col_all = repmat(weight_col, [block_num_row, 1]);

weight_block = weight_row_all+weight_col_all;

weight = kron(weight_block, ones(row_block,col_block));
end