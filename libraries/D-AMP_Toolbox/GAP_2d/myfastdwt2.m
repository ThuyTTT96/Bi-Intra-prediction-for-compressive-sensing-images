 function w = myfastdwt2(theta_temp,row,col)
 
 level = 3;
 qmf   = MakeONFilter('Daubechies',8); 
 sig_level_row = log2(row); 
 sig_level_col = log2(col); 
 T_row = get_waveletMatrix(qmf,sig_level_row,level,level);
 T_col = get_waveletMatrix(qmf,sig_level_col,level,level);

 w = reshape( shiftdim(T_col*shiftdim(T_row*reshape(theta_temp,[row col]),1),1),  [row*col 1]);

 end
