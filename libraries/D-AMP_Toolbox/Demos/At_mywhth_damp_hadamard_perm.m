function y = At_mywhth_damp_hadamard_perm(z, n,P)

 z_len = length(z);
 y = fwht([z(:); zeros(n-z_len,1)]);
 y = y(P);
end