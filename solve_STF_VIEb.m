function [u, uext] = solve_STF_VIEb(m_Gamma, m_Omega, pts, in, k0, k, ikO, kx, ax, grad_ax, tol, theta)



Omega = dom(m_Omega, 3);
Gamma = dom(m_Gamma, 3);


Uh = fem(m_Gamma, 'P1');
I = integral(Gamma, Uh, Uh);
one = ones(1, length(m_Gamma.vtx));

aO = one * I * ax(m_Gamma.vtx) / sum(ndv(m_Gamma));

kO = one * I * kx(m_Gamma.vtx) / sum(ndv(m_Gamma));
kO = (kO + 1i*ikO) / sqrt(aO);

beta = @(X) (kx(X).^2) / aO - kO^2;
alpha = @(X) (ax(X) - aO) / aO;

%%
Uh = fem(m_Omega, 'P1');


[A0, A1]     = build_STFbb(m_Gamma, k0, kO, aO, ax, grad_ax, tol);
[Av, Ad, An] = build_VIObb(m_Gamma, m_Omega, kO, beta, aO, ax, grad_ax, alpha, tol);
[S, D]       = build_potentials(m_Gamma, m_Omega, kO, tol);
F            = build_rhsb(m_Gamma, m_Omega, k0, k, theta, ax, aO); 

A = A0+A1;
traceAv = [Ad; An];
pots = [D -S];
Id = integral(Omega, Uh, ax, Uh)/aO;



M = [A traceAv;pots Id-Av];

sol = M \ F;

phi     = sol(1:size(m_Gamma, 1)); sol(1:size(m_Gamma, 1)) = [];
lambda  = sol(1:size(m_Gamma, 1)); sol(1:size(m_Gamma, 1)) = [];
u       = sol;


% uext = compute_fields(phi, lambda, u, m_Gamma, m_Omega, pts, in, k0, kO, kx, aO, ax, theta);
uext=0;
end