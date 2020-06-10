clear
clf
set (gcf, "papersize", [5,4])
set (gcf, "paperposition", [0, 0, 5, 4])
set(0,'DefaultTextInterpreter','latex')

% Constantes biologiques du système
a = 0.4; % désavantage sélectif de u, induisant la bistabilité si a>0.5 (état intermédiaire à (2a-1)/a)
K = 2; % capacité de charge de la population totale
s = .4*K; % éventuel seuil Allee de la population totale (voir lignes 40-42)

% Constantes du schéma numérique
T = 100; % temps final
L = 100; % largeur du domaine spatial
Mt = 10000; % nombre de pas temporels
Mx = 10000; % nombre de pas spatiaux
dt = T/Mt; % pas temporel
dx = L/Mx; % pas spatial
X = [0:Mx]*dx; % domaine spatial discrétisé
A = spdiags([ones(Mx+1,1) -2*ones(Mx+1,1) ones(Mx+1,1)], [-1, 0, 1], Mx+1, Mx+1); % laplacien 1D
A(1,1) = A(1,1)+1; % condition de Neumann à gauche
A(end,end) = A(end,end)+1; % condition de Neumann à droite
Aimp = eye(Mx+1)-(dt/dx^2)*A;
B = spdiags([-ones(Mx+1,1) ones(Mx+1,1)], [0, 1], Mx+1, Mx+1); % gradient 1D
B(end,end-1) = -1;
B(end,end) = 1;

% Initialisation
U = [zeros(9*Mx/20+1,1);0.99*ones(2*Mx/20,1);zeros(9*Mx/20,1)]; % u à t=0
N = [K*ones(Mx+1,1)]; % n à t=0
plot(X,U,'-b',X,N,'-r')
axis([0 Mx*dx 0 max(K,1)+0.1])
xlabel('$x$','fontsize',10)
ylabel('$(u,n)$','fontsize',10)
grid on
%print('evolution_0','-dpdflatex')

% Boucle
for i=[1:Mt]
  % Ci-dessous, choix de la dynamique de pop. sur n
  %w = (K-N).*(N-s)+1; % effet Allee fort sur n, drive d’éradication ssi (K-s)^2<4a/(1-a)
  w = K-N+1; % pas d’effet Allee sur n, drive d’éradication ssi K<a/(1-a)
  % Puis itération de la boucle
	u = Aimp\(U + 2*(dt/dx^2)*(B*log(N)).*(B*U) + dt*a*w.*U.*(1-U).*(U-(2*a-1)/a)); 
	n = Aimp\(N + dt*N.*(w.*(1+a*U.^2-2*a*U)-1)); 
	U = u;
	N = n;
	if mod(i,Mt/1000)==0
		plot(X,U,'-b',X,N,'-r')
    axis([0 Mx*dx 0 max(K,1)+0.1])
    xlabel('$x$','fontsize',10)
    ylabel('$(u,n)$','fontsize',10)
    grid on
		drawnow;
	endif  
  %if mod(i,Mt/10)==0
    %print(['evolution_',num2str(i*10/Mt)],'-dpdflatex')
  %endif
end
