function F = build_rhs_unique(mesh, k0, ax, phi, tol)


Gamma = dom(mesh, 3);
Vh = fem(mesh, 'P0');
Uh = fem(mesh, 'P1');

alpha = @(X) ax(X) - 1;


Gxy = @(X,Y) femGreenKernel(X, Y, '[H0(kr)]', k0);
dyGxy{1} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]1', k0);
dyGxy{2} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]2', k0);
dyGxy{3} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]3', k0);

dxGxy{1} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]1', k0);
dxGxy{2} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]2', k0);
dxGxy{3} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]3', k0);


I = integral(Gamma, Vh, Uh);

K = (1i/4) * integral(Gamma, Gamma, Vh, dyGxy, ntimes(Uh));
K = K - 1/(2*pi) * regularize(Gamma, Gamma, Vh,'grady[log(r)]',ntimes(Uh));

W = 1i/4 * (k0^2 * integral(Gamma, Gamma, ntimes(Uh),Gxy,ntimes(Uh)) ...
    - integral(Gamma, Gamma, nxgrad(Uh), Gxy, nxgrad(Uh)));
W = W -1/(2*pi) * (k0^2 * regularize(Gamma, Gamma, ntimes(Uh),'[log(r)]',ntimes(Uh)) ...
    - regularize(Gamma, Gamma, nxgrad(Uh),'[log(r)]',nxgrad(Uh)));

dim = size(mesh.vtx, 1);
Alpha = spdiags(alpha(mesh.vtx), 0, dim, dim);

neu = W * Alpha * phi;
dir = (0.5*I - K) * Alpha * phi;

F = [neu; dir];

end