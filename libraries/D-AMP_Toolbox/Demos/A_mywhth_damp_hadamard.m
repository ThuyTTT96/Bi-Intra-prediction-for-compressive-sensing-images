function y = A_mywhth_damp_hadamard(z, m)

 temp = fwht(z);
 y = temp(1:m);

end