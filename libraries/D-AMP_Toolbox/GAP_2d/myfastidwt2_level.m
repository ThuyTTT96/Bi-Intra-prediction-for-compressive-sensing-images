 function theta = myfastidwt2_level(w_temp,row,col,level)
 if nargin<4
 level = 2;
 end
 qmf   = MakeONFilter('Daubechies',4); 
 sig_level_row = log2(row); 
 sig_level_col = log2(col); 
 T_row = get_waveletMatrix(qmf,sig_level_row,level,level);
 T_col = get_waveletMatrix(qmf,sig_level_col,level,level);

 theta = reshape( shiftdim(T_col'*shiftdim(T_row'*reshape(w_temp,[row col]),1),1),  [row*col 1]);

 end
