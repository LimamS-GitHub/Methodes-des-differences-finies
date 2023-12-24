clear all
close all

PLOT= 1;
FREQ=10;

% 1.Initialisation de la structure EDP
EDP.a=0;EDP.b=2*pi;
EDP.t0=0;EDP.T=2;
EDP.nu=2;
k=5;
EDP.uex=@(t,x) cos(k*t)*cos(x);
EDP.f=@(t,x) -k*sin(k*t)*cos(x)+EDP.nu*cos(k*t)*cos(x);
EDP.u0=@(x) EDP.uex(EDP.t0,x);
EDP.ua=@(t) EDP.uex(t,EDP.a);
EDP.ub=@(t) EDP.uex(t,EDP.b);


% 2. Parametres de discretisation
Nx=500;
Nt=5000;

% 3. Resolution par le schema d'EULER explicite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Fonction non fournie : A IMPLEMENTER
disp('Calcul en cours...')
tic
[t,x,u]=EulerImplicite_sparse(EDP,Nt,Nx);
temps_creuse = toc
##tic
##[t,x,u]=EulerImplicite(EDP,Nt,Nx);
##temps = toc
disp('Fin du calcul.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Uex=CalculF(EDP.uex,t,x); % Solution exacte
MIN=min(min(Uex));
MAX=max(max(Uex));
% 4. Representations graphiques
if PLOT
    figure(1)
    PlotSol(t,x,u,'freq',FREQ,'title','Sol. Appr.','axis',[x(1) x(end) MIN MAX],'pause',0.1)

end

Err=abs(u-Uex);
MAX=max(max(Err));
if MAX >1
    AXIS=[];,
else
    AXIS=[x(1) x(end) 0 MAX];
end
if PLOT
    figure(2)
    PlotSol(t,x,Err,'freq',FREQ,'title','Erreur','axis',[x(1) x(end) 0 MAX],'pause',0.1)
end

Ninf=NormInf(Err);
if PLOT
  figure(3)
  plot(t,Ninf);
  xlabel('t')
  title('Erreur en norme L^\infty en espace')
end

s=sprintf('Erreur relative  (max en temps et espace) : %e',max(Ninf)/max(max(abs(Uex))));
disp(s)


