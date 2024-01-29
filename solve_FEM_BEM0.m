function u = solve_FEM_BEM0(m_Omega, ax, theta, f)



Theta = @(X) mod(atan2(X(:, 2), X(:, 1)), 2*pi);
% Theta = @(X) atan2(X(:, 2), X(:, 1));
Rad = @(X) X(:, 1).^2 + X(:, 2).^2;

dRad{1} = @(X) 2 * X(:, 1);
dRad{2} = @(X) 2 * X(:, 2);
dRad{3} = @(X) 0 * X(:, 3);

dTheta{1} = @(X) -sin(Theta(X)) ./ Rad(X).^(1/2);
dTheta{2} = @(X)  cos(Theta(X)) ./ Rad(X).^(1/2);
dTheta{3} = @(X) 0*X(:, 1);



Uin = @(X) Rad(X).^(1/3) .* sin(2 / 3 * Theta(X));
Uex = @(X) 0.5 * log( (X(:, 1) + 1/4).^2 + (X(:, 2) + 1/4).^2);


gradxUin{1} = @(X) sin(2 / 3 * Theta(X)) .* dRad{1}(X) .* Rad(X).^(-2/3) /3 + ...
                   Rad(X).^(1/3) .* cos(2 / 3 * Theta(X)) .* dTheta{1}(X) * 2 / 3;
               
gradxUin{2} = @(X) sin(2 / 3 * Theta(X)) .* dRad{2}(X) .* Rad(X).^(-2/3) /3 + ...
                   Rad(X).^(1/3) .* cos(2 / 3 * Theta(X)) .* dTheta{2}(X) * 2 / 3;
gradxUin{3} = @(X) 0 .* X(:, 1);


gradxUex{1} = @(X) (X(:, 1) + 1/4) ./ ((X(:, 1) + 1/4).^2 + (X(:, 2) + 1/4).^2);
gradxUex{2} = @(X) (X(:, 2) + 1/4) ./ ((X(:, 1) + 1/4).^2 + (X(:, 2) + 1/4).^2);
gradxUex{3} = @(X) 0 .* X(:, 1);


uinc = @(X) Uex(X) - Uin(X);

gradxUinc{1} = @(X) gradxUex{1}(X) - gradxUin{1}(X);
gradxUinc{2} = @(X) gradxUex{2}(X) - gradxUin{2}(X);
gradxUinc{3} = @(X) 0 .* X(:, 1);

k0 = 1;
% uinc = @(X) exp(1i*k0*(cos(theta)*X(:,1) +sin(theta)*X(:,2)));
% X0 = [0.2 0.2];
% uinc = @(X) 1i / 4 * besselh(0, k0 * sqrt( (X(:, 1) - X0(1)).^2 + (X(:, 2) - X0(2)).^2) );
% uinc = @(X) 1 / 2 / pi * log(sqrt( (X(:, 1) - X0(1)).^2 + (X(:, 2) - X0(2)).^2));
% uinc = @(X) ones(size(X, 1), 1);
m_Gamma = bnd(m_Omega);


Omega = dom(m_Omega, 3);
Gamma = dom(m_Gamma, 3);

%%


Gxy = @(X,Y) femGreenKernel(X, Y, '[log(r)]', 0);
dyGxy{1} = @(X,Y) femGreenKernel(X, Y, 'grady[log(r)]1', 0);
dyGxy{2} = @(X,Y) femGreenKernel(X, Y, 'grady[log(r)]2', 0);
dyGxy{3} = @(X,Y) femGreenKernel(X, Y, 'grady[log(r)]3', 0);

%% Boundary Integral Operators
Vh = fem(m_Gamma, 'P0');
Uh = fem(m_Gamma, 'P1');


V = -1/2/pi * integral(Gamma, Gamma, Vh, Gxy, Vh);
V = V - 1/2/pi * regularize(Gamma, Gamma, Vh,'[log(r)]',Vh);

K = -1/2/pi * integral(Gamma, Gamma, Vh, dyGxy, ntimes(Uh));
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


Uh = fem(m_Omega, 'P1');


[~,I1,I2] = intersect(Wh.unk,Uh.unk,'rows','stable');
P = sparse(I1,I2,1,length(Wh),length(Uh));



A = A;

B = -integral(Gamma, Uh, Vh);

C = (0.5* Id + K) * P;

D = V;

T = [A B; C D];


Neu = M0 \ Neu;
Dir = M1 \ Dir;

Rhs = (0.5 * Id  + K) * Dir;

F = integral(Omega, Uh, f) - integral(Gamma, Uh, Vh) * Neu;
Rhs = [F; Rhs];

% 
% Rhs = integral(Gamma, Vh, uinc);
% Rhs = [zeros(size(A, 1), 1); Rhs];
% 
sol = T \ Rhs;
u      = sol(1:size(A, 1)); sol(1:size(A,1)) = [];
lambda = sol; 

end