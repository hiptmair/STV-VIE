function [S, D] = build_potentials2(mesh, pts, k0)

Uh = fem(mesh, 'P1');
Gamma = dom(mesh, 3);
Omega = dom(pts, 3);

[Xq, Wq, elt2qud] = Omega.qud; 

Wq = spdiags(Wq, 0, size(Wq, 1), size(Wq, 1));

% in = inpolygon(pts.vtx(:, 1), pts.vtx(:, 2), mesh.vtx(:, 1), mesh.vtx(:, 2));

Gxy = @(X,Y) femGreenKernel(X, Y, '[H0(kr)]', k0);
dyGxy{1} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]1', k0);
dyGxy{2} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]2', k0);
dyGxy{3} = @(X,Y) 0*X(:, 1);

dxGxy{1} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]1', k0);
dxGxy{2} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]2', k0);
dxGxy{3} = @(X,Y) 0*X(:, 1);



% disp('Assembly of Layer Potentials');

% Single Layer

S0   = 1i/4 * integral(Xq, Gamma, Gxy, Uh);
Sr0  = -1/(2*pi) .* regularize(Xq, Gamma,'[log(r)]', Uh);
S   = S0 + Sr0;


% Double layer

D0  = 1i/4 * integral(Xq, Gamma, dyGxy, ntimes(Uh));
Dr0 = -1/(2*pi) .* regularize(Xq, Gamma,'grady[log(r)]', ntimes(Uh));
D  = D0 + Dr0;


Omega = dom(pts, 3);
Uh = fem(pts, 'P1');


M = femLagrangePn(Uh, Omega); 
[~,P] = Uh.unk;    
M = M * P;  

S = M' * Wq * S;
D = M' * Wq * D;


end