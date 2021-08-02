function y = A_mywhth_damp_hadamard_perm(z, m, Q)

 temp = fwht(z(Q));
 y = temp(1:m);

end