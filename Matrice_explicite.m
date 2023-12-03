function M = Matrice_explicite (d,nu,delta_t,delta_x)
  M = eye(d)*(-2+delta_x^2/(delta_t*nu));
  for i =1: d-1
    M(i,i+1) = 1;
    M(i+1,i) = 1;
  endfor
end

