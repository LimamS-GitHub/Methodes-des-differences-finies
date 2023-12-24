function K=spMatTriDiagSca (v,u,w)
  d = length(u);
  K=sparse (d,d) ;
  for i=1:d-1
    K (i,i)= u(i);
    K (i+1,i)= w(i);
    K (i,i+1)= v(i);
  end
  K(d,d)=u(d);
end

