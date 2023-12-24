function [t,x,u]=EulerImplicite_sparse(EDP,Nt,Nx)
  delta_t = (EDP.T-EDP.t0)/Nt;
  delta_x = (EDP.b-EDP.a)/Nx;
  t = EDP.t0:delta_t:EDP.T+EDP.t0;
  x = EDP.a:delta_x:EDP.b;
  u = zeros (Nx+1,Nt+1);

  u(:,1) = EDP.u0(x);
  i = ones(1,Nx-1);
  k = EDP.nu*delta_t/delta_x^2;
  w = ones(1,Nx-1)*(1+2*k);
  v = ones(1,Nx-2)*(-k);

  Ah = spMatTriDiagSca (v,w,v)
  A = sparse(Nx+1,Nx+1);
  A(2:Nx,2:Nx) = Ah;
  A(1,1)=EDP.delta+3*EDP.mua/(2*delta_x);
  A(1,2)=-2*EDP.mua/delta_x;
  A(1,3)=EDP.mua/(2*delta_x);
  A(end,end)=EDP.deltb+3*EDP.mub/(2*delta_x);
  A(end,end-1)= -2*EDP.mub/delta_x;
  A(end,end-2)= EDP.mub/(2*delta_x);
  A(2,1) = -k;
  A(end-1,end) = -k;
  F=delta_t*CalculF(EDP.f,t(1:Nt+1),x(2:Nx));

  f = sparse (Nx+1,Nt+1);
  f(2:Nx,1:Nt+1)=F;
  f(1,:)=EDP.ua(t(1:Nt+1));
  f(end,:)=EDP.ub(t(1:Nt+1));
  for n = 1:Nt
    u(1:Nx+1,n+1) = (A)\(f(:,n+1)+[0;u(2:Nx,n);0]);
  end
end

