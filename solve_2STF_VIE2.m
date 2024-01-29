function [u, uext] = solve_2STF_VIE2(m_Gamma, m_Omega, pts, in, k0, k, kx, ax, tol, theta)



Omega = dom(m_Omega, 3);
Gamma = dom(m_Gamma, 3);

% disp(["Degrees of Freedom Volume: ", num2str(size(m_Omega.vtx,1))]);

% disp(["Degrees of Freedom Boundary: ", num2str(size(m_Gamma.vtx,1))]);

Uh = fem(m_Gamma, 'P1');
I = integral(Gamma, Uh, Uh);
one = ones(1, length(m_Gamma.vtx));




aO = one * I * ax(m_Gamma.vtx) / sum(ndv(m_Gamma));

kO = one * I * kx(m_Gamma.vtx) / sum(ndv(m_Gamma));
kO = kO / sqrt(aO);

beta = @(X) (kx(X).^2) / aO - kO^2;
alpha = @(X) (ax(X) - aO) / aO;


I0 = integral(Gamma, Uh, Uh);
I1 = integral(Gamma, Uh, Uh);

Id = [I0 sparse(size(I0, 1), size(I1, 2)); ...
      sparse(size(I1, 1), size(I0, 2)) I1];

%%
Uh = fem(m_Omega, 'P1');


[A0, A1]     = build_2STF2(m_Gamma, k0, kO, ax, alpha, tol);
[Av, Ad, An] = build_2VIO2(m_Gamma, m_Omega, kO, beta, ax, alpha, tol);
[S, D]       = build_potentials2(m_Gamma, m_Omega, kO);
F            = build_rhs2(m_Gamma, m_Omega, k0, k, theta, ax);

F = -F;


A = Id + A0 - A1;
traceAv = -[Ad; An];
pots = [D -S];
Id = integral(Omega, Uh, Uh);



M = [A traceAv;pots Id-Av];

% tic

% disp('--------------------------------')
% disp("Solve Linear System")

sol = M \ F;
% toc

phi     = sol(1:size(m_Gamma, 1)); sol(1:size(m_Gamma, 1)) = [];
lambda  = sol(1:size(m_Gamma, 1)); sol(1:size(m_Gamma, 1)) = [];
u       = sol;

% disp('--------------------------------')


uext = compute_fields(phi, lambda, u, m_Gamma, m_Omega, pts, in, k0, kO, kx, aO, ax, theta);

end