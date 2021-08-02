 function w = myfastdct2(theta_temp,row,col)
 base_row = dct(eye(row));
 base_col = dct(eye(col));
 w = reshape( shiftdim(base_col*shiftdim(base_row*reshape(theta_temp,[row col]),1),1),  [row*col 1]);

 end
