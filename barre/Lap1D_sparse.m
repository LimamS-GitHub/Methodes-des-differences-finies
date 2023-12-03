function K=Lap1D_sparse ( d , x , nu )
  K=sparse (d , d ) ;
  K (1,1) =-2*nu(x(1));K (1,2) = nu(x(1));
  for i=2:d-1
    K (i,i)=-2*nu(x(i));
    K (i,i-1)=nu(x(i));
    K (i,i+1)=nu(x(i));
  end
  K (d,d)=-2*nu(x(d));
  K (d,d-1)=nu(x(d-1));
end
