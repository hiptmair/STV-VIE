function [u, uext] = solve_STF_VIE(m_Gamma, m_Omega, pts, in, k0, kx, tol, theta)

tic

Omega = dom(m_Omega, 3);
Gamma = dom(m_Gamma, 3);

% disp(["Degrees of Freedom Volume: ", num2str(size(m_Omega.vtx,1))]);

% disp(["Degrees of Freedom Boundary: ", num2str(size(m_Gamma.vtx,1))]);

Uh = fem(m_Gamma, 'P1');
I = integral(Gamma, Uh, Uh);
one = ones(1, length(m_Gamma.vtx));

kO = one * I * kx(m_Gamma.vtx) / sum(ndv(m_Gamma));

% kO = one * I / sum(ndv(m_Gamma));
% kO = 6;
beta = @(X) kx(X).^2 - kO^2;

%%
Uh = fem(m_Omega, 'P1');


A            = build_STF(m_Gamma, k0, kO, tol);toc
[Av, traceA] = build_VIO(m_Gamma, m_Omega, k0, kO, beta, tol);toc
pots         = build_potentials(m_Gamma, m_Omega, kO, tol);toc
[F, Fg]      = build_rhs(m_Gamma, m_Omega, k0, theta); toc


% A = A0+A1;
% traceAv = [Ad; An];
% traceAv = [An; Ad];


% pots = [D -S];
Id = integral(Omega, Uh, Uh);


% 
% M = bmm({A0.blk{1,1} + A1.blk{1,1}, A0.blk{1,2} + A1.blk{1,2}, An; ...
%          A0.blk{2,1} + A1.blk{2,1}, A0.blk{2,2} + A1.blk{2,2}, Ad; ...
%          D, -S, Id-Av});
%           
% M = [(A0.blk{1,1}+A1.blk{1,1}) (A0.blk{1,2} + A1.blk{1,2}) An; ...
%          (A0.blk{2,1} + A1.blk{2,1}) (A0.blk{2,2}+A1.blk{2,2}) Ad; ...
%          D -S Id-Av];
     

M = [A traceA; ...
     pots Id-Av];

% M = [A traceAv;pots Id-Av];
% M = bmm({A, traceAv; pots, Id-Av});

% M0 = bmm({A, sparse(0*traceAv);sparse(0*pots), Id});
% 
% A0 = [A0.blk{1,1} A0.blk{1,2};A0.blk{2,1} A0.blk{2,2}];
% A1 = [A1.blk{1,1} A1.blk{1,2};A1.blk{2,1} A1.blk{2,2}];
    

 
% M0 = bmm({A, sparse(size(Ad, 1)+size(An,1), size(Ad, 2)); ...
%           sparse(size(S, 1), size(S,2)+ size(D, 2)), Id});

 % tic

% disp('--------------------------------')
% disp("Solve Linear System")


toc
disp(" ");
disp("Direct Solver")
tic
sol = M \ F;
toc

disp(" ");
disp("Iterative Solver: no preconditioner")
tic
sol2 = gmres(@(x) M*x, F, [], 1e-6, 500);
toc

disp(" ");
disp("Iterative Solver: precond. M0 inverse")
tic

% [L0, U0] = lu(A);
M0 = bmm({A, zeros(traceA); ...
      zeros(pots), speye(size(Id))});

[L0, U0] = lu(M0);

toc
sol3 = gmres(@(x) M*x, F, [], 1e-6, 500, @(x) L0 \ x, @(x) U0 \ x);
toc


phi     = sol(1:size(m_Gamma, 1)); sol(1:size(m_Gamma, 1)) = [];
lambda  = sol(1:size(m_Gamma, 1)); sol(1:size(m_Gamma, 1)) = [];
u       = sol;

% disp('--------------------------------')

uext =0;
% uext = compute_fields(phi, lambda, u, m_Gamma, m_Omega, pts, in, k0, kO, kx, theta);

end