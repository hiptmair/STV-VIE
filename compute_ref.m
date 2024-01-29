function u = compute_ref(mref, m, uh, k0, kx, theta, tol, ax)

if nargin < 8
   ax = @(X) ones(size(X, 1), 1); 
end

Omega = dom(mref, 3);
Vh = fem(mref, 'P1');


G = @(X, Y) femGreenKernel(X, Y, '[H0(kr)]', k0);

dxGxy{1} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]1', k0);
dxGxy{2} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]2', k0);
dxGxy{3} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]3', k0);


comm_dxGxy{1} = @(X,Y) (ax(Y) - ax(X)) .* femGreenKernel(X, Y, 'gradx[H0(kr)]1', k0);
comm_dxGxy{2} = @(X,Y) (ax(Y) - ax(X)) .* femGreenKernel(X, Y, 'gradx[H0(kr)]2', k0);
comm_dxGxy{3} = @(X,Y) (ax(Y) - ax(X)) .* femGreenKernel(X, Y, 'gradx[H0(kr)]3', k0);





uinc = @(X) exp(1i*k0*(cos(theta)*X(:,1) +sin(theta)*X(:,2)));
beta = @(X) kx(X).^2 - k0^2;
alpha = @(X) ax(X) - 1;

dim = size(mref.vtx, 1);
Beta = spdiags(beta(mref.vtx), 0, dim, dim); 
% tic
K = (1i/4) * integral(m.vtx, Omega, G, Vh, tol);
Kreg = -(1/(2 * pi)) * regularize(m.vtx, Omega, '[log(r)]',Vh);

% toc

aux = uinc(m.vtx);


K1 = (1i / 4) * integral(m.vtx, Omega, dxGxy, grad(Vh));
commK = (1i / 4) * integral(m.vtx, Omega, comm_dxGxy, grad(Vh));


reg = - 1 / (2 * pi) * regularize(m.vtx, Omega, 'grady[log(r)]', grad(Vh));

dim = size(m.vtx, 1);
Alpha = spdiags(alpha(m.vtx),0,dim,dim);


u = aux + (K + Kreg)* (Beta*uh) + (Alpha * (K1 - reg) + commK) * uh;

end