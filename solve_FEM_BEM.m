function [u, lambda] = solve_FEM_BEM(m_Omega, k0, kx, theta, k, ax, f)
tic
if nargin < 5
   k = 1;
   ax = @(X) ones(size(X, 1), 1);
   f = @(X) zeros(size(X, 1), 1);
end

k2 = @(X) kx(X).^2;

[uinc, gradxUinc] = incident_field(k0, k, ax, theta, 'PW');


m_Gamma = bnd(m_Omega);


Omega = dom(m_Omega, 3);
Gamma = dom(m_Gamma, 3);

%%


Gxy = @(X,Y) femGreenKernel(X, Y, '[H0(kr)]', k0);
dyGxy{1} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]1', k0);
dyGxy{2} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]2', k0);
dyGxy{3} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]3', k0);

%% Boundary Integral Operators
Vh = fem(m_Gamma, 'P0');
Uh = fem(m_Gamma, 'P1');


V = (1i/4) * integral(Gamma, Gamma, Vh, Gxy, Vh);
V = V - 1/2/pi * regularize(Gamma, Gamma, Vh,'[log(r)]',Vh);

K = (1i/4) * integral(Gamma, Gamma, Vh, dyGxy, ntimes(Uh));
K = K - 1/(2*pi) * regularize(Gamma, Gamma, Vh,'grady[log(r)]',ntimes(Uh));

Wh = fem(m_Gamma,'P1');

Id = integral(Gamma, Vh, Uh);

M0 = integral(Gamma, Vh, Vh);
M1 = integral(Gamma, Uh, Uh);

Dir = integral(Gamma, Uh, uinc);
Neu = integral(Gamma, ntimes(Vh), gradxUinc);

%% Finite Element Matrices
Uh = fem(m_Omega, 'P1');

A = integral(Omega,grad(Uh), ax, grad(Uh));
M = integral(Omega, Uh, k2, Uh);

% P= restriction(Uh,m_Gamma); % operator restriction on Gamma


Uh = fem(m_Omega, 'P1');


[~,I1,I2] = intersect(Wh.unk,Uh.unk,'rows','stable');
P = sparse(I1,I2,1,length(Wh),length(Uh));



A = A - M;

B = -integral(Gamma, Uh, Vh);

C = (0.5* Id + K) * P;

D = V;

T = [A B; C D];


Neu = M0 \ Neu;
Dir = M1 \ Dir;

Rhs = (0.5 * Id  + K) * Dir;

F = integral(Omega, Uh, f) - integral(Gamma, Uh, Vh) * Neu;

% Rhs = Id * Dir;
% F = integral(Omega, Uh, f);

Rhs = [F; Rhs];



toc
disp(" ");
disp("Direct Solver")
tic

sol = T \ Rhs;
toc


% disp(" ");
% disp("Iterative Solver: no preconditioner")
% tic
% sol2 = gmres(T, Rhs, [], 1e-6, 500);
% toc




u      = sol(1:size(A, 1)); sol(1:size(A,1)) = [];
lambda = sol; 

end