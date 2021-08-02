function y = At_mywhth_damp_hadamard(z, n)

 z_len = length(z);
 y = fwht([z; zeros(n-z_len,1)]);

end