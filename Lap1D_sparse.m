function K=Lap1D_sparse ( d )
  K=sparse (d , d ) ;
  K (1,1) =-2;K (1,2) = 1;
  for i=2:d-1
    K (i,i)=-2;
    K (i,i-1)=1;
    K (i,i+1)=1;
  end
  K (d,d)=-2;
  K (d,d-1)=1;
end
