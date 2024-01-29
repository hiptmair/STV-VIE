function [F, Fg] = build_rhs(mesh,pts, k0, theta, ax)

if nargin < 5
   ax = @(X) ones(size(X, 1), 1); 
end

Gamma = dom(mesh, 3);

Vh = fem(mesh, 'P0');
Uh = fem(mesh, 'P1');

uinc = @(X) exp(1i*k0*(cos(theta)*X(:,1) +sin(theta)*X(:,2)));
gradxUinc{1} = @(X) (1 ./ ax(X)) .* 1i*k0*cos(theta) .* uinc(X);
gradxUinc{2} = @(X) (1 ./ ax(X)) .* 1i*k0*sin(theta) .* uinc(X);
gradxUinc{3} = @(X) (1 ./ ax(X)) .* 1i*k0*0 .* uinc(X);

dir = integral(Gamma, Vh, uinc);
neu = integral(Gamma, ntimes(Uh), gradxUinc);

% F = [dir; neu; zeros(length(pts.vtx), 1)];


F = [neu; dir; zeros(length(pts.vtx), 1)];

Fg = [neu; dir];
end