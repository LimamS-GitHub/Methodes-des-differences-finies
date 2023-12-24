function M = Lap1D (d)
  M = -2*eye(d);
  for i = 1:d-1
    M (i,i+1) = 1;
    M (i+1,i) = 1;
  end
end
