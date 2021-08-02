 function theta = myfastidct2(w_temp,row,col)
 base_row = dct(eye(row));
 base_col = dct(eye(col));
 theta = reshape( shiftdim(base_col'*shiftdim(base_row'*reshape(w_temp,[row col]),1),1),  [row*col 1]);

 end
