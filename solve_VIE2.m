function [uh, uext] = solve_VIE2(mesh, pts, k0, k, kx, ax, typ, tol, theta)

Omega = dom(mesh, 3);

m_Omega = mesh;
m_Gamma = bnd(mesh);

Gamma = dom(m_Gamma, 3);

G = @(X, Y) femGreenKernel(X, Y, '[H0(kr)]', k0);

dxGxy{1} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]1', k0);
dxGxy{2} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]2', k0);
dxGxy{3} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]3', k0);


comm_dxGxy{1} = @(X,Y) (ax(Y) - ax(X)) .* femGreenKernel(X, Y, 'gradx[H0(kr)]1', k0);
comm_dxGxy{2} = @(X,Y) (ax(Y) - ax(X)) .* femGreenKernel(X, Y, 'gradx[H0(kr)]2', k0);
comm_dxGxy{3} = @(X,Y) (ax(Y) - ax(X)) .* femGreenKernel(X, Y, 'gradx[H0(kr)]3', k0);
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uh = fem(m_Gamma, 'P1');
% Vh = fem(m_Gamma, 'P0');
% 
% M0 = integral(Gamma, Vh, Vh);
% M1 = integral(Gamma, Uh, Uh);
% 
% Dir = integral(Gamma, Uh, uinc);
% Neu = integral(Gamma, ntimes(Vh), gradxUinc);
% 
% Neu = M0 \ Neu;
% Dir = M1 \ Dir;
% 
% [S, D]       = build_potentials(m_Gamma, m_Omega, k0);
% 
% F = -S * Neu + D * Dir;
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beta = @(X) kx(X).^2 - k0^2;
alpha = @(X) ax(X) - 1.0;


[Xq, Wq, elt2qud] = Omega.qud;
Wq = spdiags(Wq, 0, size(Wq, 1), size(Wq, 1));

Vh = fem(mesh, typ);

uinc = @(X) exp(1i*k0*(cos(theta)*X(:,1) +sin(theta)*X(:,2)));
F = integral(Omega, Vh, uinc);

K = (1i/4) * integral(Omega, Omega,  Vh, G, Vh, tol);
Kreg = -(1/(2 * pi)) * regularize(Omega,Omega,Vh,'[log(r)]',Vh);

dim = size(m_Omega.vtx, 1);
Beta = spdiags(beta(m_Omega.vtx), 0, dim, dim); 
Betah = hmx(double(m_Omega.vtx),double(m_Omega.vtx),Beta,tol);

% Av = (K + Kreg) * Beta;
Avh = (K + Kreg) * Betah;

K = (1i / 4) * integral(Xq, Omega, dxGxy, grad(Vh), tol);
commK = (1i / 4) * integral(Xq, Omega, comm_dxGxy, grad(Vh), tol);


reg = - 1 / (2 * pi) * regularize(Xq, Omega, 'grady[log(r)]', grad(Vh));

m = size(Xq, 1);
Alpha = spdiags(alpha(Xq),0,m,m);
Alphah = hmx(double(Xq),double(Xq),Alpha,tol);

Mv = femLagrangePn(Vh, Omega); 

[~,Pv] = Vh.unk;  

Mv = Mv * Pv;


% Av = Av + Mv' * Wq * (Alpha * (K - reg) + commK);


Avh = Avh + Mv' * Wq * (Alphah * (K - reg) + commK);
I = integral(Omega, Vh, Vh);
% M = I - Av;

Mh = I - Avh;


% sol = M \ F;
solh = Mh \ F;

% uh = sol;
uh = solh;


% Evaluates using Representation formula

K = (1i/4) * integral(pts.vtx, Omega, G, Vh);
Kreg = -(1/(2 * pi)) * regularize(pts.vtx, Omega, '[log(r)]',Vh);

aux = uinc(pts.vtx);


uext = aux + (K + Kreg)* (Beta*uh);


K = (1i / 4) * integral(pts.vtx, Omega, dxGxy, grad(Vh), tol);
commK = (1i / 4) * integral(pts.vtx, Omega, comm_dxGxy, grad(Vh), tol);


reg = - 1 / (2 * pi) * regularize(pts.vtx, Omega, 'grady[log(r)]', grad(Vh));

m = size(pts.vtx, 1);
Alpha = spdiags(alpha(pts.vtx),0,m,m);



uext = uext + (Alpha * (K - reg) + commK) * uh;





end
 
 
 
 
 
 
 
 