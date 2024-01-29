function [S, D] = build_potentials(mesh, pts, k0, tol)

Vh = fem(mesh, 'P0');
Uh = fem(mesh, 'P1');
Gamma = dom(mesh, 3);
Omega = dom(pts, 3);

[Xq, Wq, elt2qud] = Omega.qud; 

Wq = spdiags(Wq, 0, size(Wq, 1), size(Wq, 1));

% in = inpolygon(pts.vtx(:, 1), pts.vtx(:, 2), mesh.vtx(:, 1), mesh.vtx(:, 2));

Gxy = @(X,Y) femGreenKernel(X, Y, '[H0(kr)]', k0);
dyGxy{1} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]1', k0);
dyGxy{2} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]2', k0);
dyGxy{3} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]3', k0);

dxGxy{1} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]1', k0);
dxGxy{2} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]2', k0);
dxGxy{3} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]3', k0);



% disp('Assembly of Layer Potentials');

Omega = dom(pts, 3);
uh = fem(pts, 'P1');
% Single Layer
% 
% S0   = 1i/4 * integral(Xq, Gamma, Gxy, Vh, tol);
% Sr0  = -1/(2*pi) .* regularize(Xq, Gamma,'[log(r)]', Vh);
% S   = S0 + Sr0;


SS0   = 1i/4 * integral(Omega, Gamma, uh, Gxy, Vh);
SSr0  = -1/(2*pi) .* regularize(Omega, Gamma,uh,'[log(r)]', Vh);
SS   = SS0 + SSr0;

% Double layer

% D0  = 1i/4 * integral(Xq, Gamma, dyGxy, ntimes(Uh), tol);
% Dr0 = -1/(2*pi) .* regularize(Xq, Gamma,'grady[log(r)]', ntimes(Uh));
% D  = D0 + Dr0;


DD0  = 1i/4 * integral(Omega, Gamma, uh, dyGxy, ntimes(Uh));
DDr0 = -1/(2*pi) .* regularize(Omega, Gamma, uh, 'grady[log(r)]', ntimes(Uh));
DD  = DD0 + DDr0;



% M = femLagrangePn(uh, Omega); 
% [~,P] = uh.unk;    
% M = M * P;  
% 
% S = M' * Wq * S;
% D = M' * Wq * D;

S = SS;
D = DD;
% norm(S-full(SS))
% norm(D-full(DD))

pots = [D -S];


end