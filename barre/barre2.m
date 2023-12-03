close all;
clear all;
L = 6;
EDP.a=0;EDP.b=L;
EDP.t0=0;EDP.T=10;
EDP.nu=@(x)1*(x<=L/3)+0.1*(x>L/3)*(x<2*L/3)+1*(x>=2*L/3);
EDP.f=@(t,x) 0;
EDP.u0=@(x) 100;
EDP.ua=@(t) ((t<=1).*(100-90*t)+(t>1)*10);
EDP.ub=@(t) (t<=1).*(100-80*t)+(t>1)*20;
Nx = Nt = 100;
[t,x,u]=EulerImplicite_sparse(EDP,Nt,Nx);
[T,X]= meshgrid(t,x);
figure;
MIN=min(min(u));
MAX=max(max(u));

if 1
    figure(1)
    PlotSol(t,x,u,'freq',1,'title','Sol. Appr.','axis',[x(1) x(end) MIN MAX],'pause',0.1)
end
figure;
surf(T,X,u);
