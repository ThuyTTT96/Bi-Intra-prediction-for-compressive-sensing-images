function weight = my_get_weight_wavelettree(row, col,block)
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
%block_num_col = length(col_block);

weight = zeros(row,col);
weight_row = 1e-5+cos((0:(block_num_row-1))/block_num_row*pi/2).^2;
%weight_col = 1e-5+cos((0:(block_num_col-1))/block_num_col*pi/2).^2;

weight((1:row_block(1)), (1:col_block(1))) = kron(weight_row(1),ones(row_block(1),col_block(1)));

rowsum = row_block(1);
colsum = col_block(1);
for nR = 2:block_num_row
    
    temp_weight =  kron(weight_row(nR), ones(row_block(nR),col_block(nR)));
    weight(rowsum+(1:row_block(nR)), colsum+(1:col_block(nR))) = temp_weight;
    weight(rowsum+(1:row_block(nR)), 1:colsum) = temp_weight;
    weight(1:rowsum, colsum+(1:col_block(nR))) = temp_weight;
    
    rowsum = rowsum + row_block(nR);
    colsum = colsum + col_block(nR);
end


end