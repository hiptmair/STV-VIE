function F = build_rhs2b(mesh,pts, k0, k, theta, ax, aO)

Gamma = dom(mesh, 3);

Uh = fem(mesh, 'P1');



[uinc, gradxUinc] = incident_field(k0, k, ax, theta, 'PW');

uinc = @(X) ax(X) ./ aO .* uinc(X);

gradxUinc{1} = @(X) (1 ./ aO) .* gradxUinc{1}(X);
gradxUinc{2} = @(X) (1 ./ aO) .* gradxUinc{2}(X);
gradxUinc{3} = @(X) (1 ./ aO) .* gradxUinc{3}(X);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dir = integral(Gamma, Uh, uinc);
neu = integral(Gamma, ntimes(Uh), gradxUinc);

F = [dir; neu; zeros(length(pts.vtx), 1)];

end