function [t,x,u]=EulerImplicite_sparse(EDP,Nt,Nx)
  delta_t = (EDP.T-EDP.t0)/Nt;
  delta_x = (EDP.b-EDP.a)/Nx;
  t = EDP.t0:delta_t:EDP.T+EDP.t0;
  x = EDP.a:delta_x:EDP.b;
  u = zeros (Nx+1,Nt+1);
  u(:,1) = EDP.u0(x);
  u(1,:) = EDP.ua(t);
  u(Nx+1,:) = EDP.ub(t);
  i = ones(1,Nx-1);
  Id = spdiags(i',0,Nx-1,Nx-1);
  Ah = Id - delta_t/delta_x^2*Lap1D_sparse(Nx-1,x,EDP.nu);
  F=CalculF(EDP.f,t(1:Nt+1),x(2:Nx));
  F(1,:) = F(1,:) + EDP.nu(x(2))/delta_x^2*EDP.ua(t(1:Nt+1));
  F(end,:) = F(end,:) + EDP.nu(x(end-1))/delta_x^2*EDP.ub(t(1:Nt+1));
  for n = 1:Nt
    u(2:Nx,n+1) = (Ah)\(u(2:Nx,n) + delta_t * F(:,n+1));
  end
end

