function gradu = compute_gradref(mref, m, uh, k0, kx, theta, tol)

Omega = dom(mref, 3);
Vh = fem(mref, 'P1');


G = @(X, Y) femGreenKernel(X, Y, '[H0(kr)]', k0);

dyGxy{1} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]1', k0);
dyGxy{2} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]2', k0);
dyGxy{3} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]3', k0);

uinc = @(X) exp(1i*k0*(cos(theta)*X(:,1) +sin(theta)*X(:,2)));
gradxUinc{1} = @(X) 1i*k0*cos(theta) .* uinc(X);
gradxUinc{2} = @(X) 1i*k0*sin(theta) .* uinc(X);
gradxUinc{3} = @(X) 1i*k0*0 .* uinc(X);

beta = @(X) kx(X).^2 - k0^2;

dim = size(mref.vtx, 1);
Beta = spdiags(beta(mref.vtx), 0, dim, dim); 
% tic
K = (1i/4) * integral(m.vtx, Omega, dyGxy{1}, Vh, tol);
Kreg = -(1/(2 * pi)) * regularize(m.vtx, Omega, 'grady[log(r)]1',Vh);

% toc

aux = gradxUinc{1}(m.vtx);

gradu{1} = aux - (K + Kreg)* (Beta*uh);


K = (1i/4) * integral(m.vtx, Omega, dyGxy{2}, Vh, tol);
Kreg = -(1/(2 * pi)) * regularize(m.vtx, Omega, 'grady[log(r)]2',Vh);

% toc

aux = gradxUinc{2}(m.vtx);

gradu{2} = aux - (K + Kreg)* (Beta*uh);
end