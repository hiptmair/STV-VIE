function [u, lambda] = solve_FEM_BEM3D(m_Omega, k0, kx, dir, ax, f, tol)
tic
if nargin < 5
   ax = @(X) ones(size(X, 1), 1);
   f = @(X) zeros(size(X, 1), 1);
   tol = 1e-6;
end

k2 = @(X) kx(X).^2;

[uinc, gradxUinc] = incident_field_3D(k0, dir, 'PW');


m_Gamma = bnd(m_Omega);


Omega = dom(m_Omega, 4);
Gamma = dom(m_Gamma, 3);

%%


Gxy = @(X,Y) femGreenKernel(X, Y, '[exp(ikr)/r]', k0);
dyGxy{1} = @(X,Y) femGreenKernel(X, Y, 'grady[exp(ikr)/r]1', k0);
dyGxy{2} = @(X,Y) femGreenKernel(X, Y, 'grady[exp(ikr)/r]2', k0);
dyGxy{3} = @(X,Y) femGreenKernel(X, Y, 'grady[exp(ikr)/r]3', k0);

%% Boundary Integral Operators
Vh = fem(m_Gamma, 'P0');
Uh = fem(m_Gamma, 'P1');


V = (1/4/pi) * integral(Gamma, Gamma, Vh, Gxy, Vh, tol);
V = V + (1/4/pi) * regularize(Gamma, Gamma, Vh,'[1/r]',Vh);

K = (1/4/pi) * integral(Gamma, Gamma, Vh, dyGxy, ntimes(Uh), tol);
K = K + (1/4/pi) * regularize(Gamma, Gamma, Vh,'grady[1/r]',ntimes(Uh));

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



Uh = fem(m_Omega, 'P1');


[~,I1,I2] = intersect(Wh.unk,Uh.unk,'rows','stable');
P = sparse(I1,I2,1,length(Wh),length(Uh));


KK = (1/4/pi) * integral(Gamma, Gamma, Vh, dyGxy, ntimes(Uh), tol);
KK = KK + (1/4/pi) * regularize(Gamma, Gamma, Vh,'grady[1/r]',ntimes(Wh)) * P;


IId = integral(Gamma, Vh, Uh);

A = A - M;

B = -integral(Gamma, Uh, Vh);

C = 0.5 * IId - KK;

D = V;
toc
% T = [A B; C D];

T = bmm({A, B; C, D});


Neu = M0 \ Neu;
Dir = M1 \ Dir;

Rhs = (0.5 * Id  - K) * Dir;

F = integral(Omega, Uh, f) + integral(Gamma, Uh, Vh) * Neu;


% Rhs = Dir;
% F = integral(Omega, Uh, f);
Rhs = [F; Rhs];

sol = T \ Rhs;




u      = sol(1:size(A, 1)); sol(1:size(A,1)) = [];
lambda = sol; 

end