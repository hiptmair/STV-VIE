function uext = compute_fields(phi, lambda, u, m_Gamma, m_Omega, pts, in, k0, kO, kx, theta, aO, ax)

if nargin < 13
   aO = 1;
   ax = @(X) ones(size(X,1), 1);
end

if sum(in) ~= size(in, 1)

beta = @(X) (kx(X).^2) / aO - kO^2;

alpha = @(X) (ax(X) - aO) / aO;
check = norm(beta(m_Omega.vtx));

check2 = norm(alpha(m_Omega.vtx));

% disp("Compute Exterior and Interior Field");

Gamma = dom(m_Gamma, 3);
Omega = dom(m_Omega, 3);


G = @(X, Y) femGreenKernel(X, Y, '[H0(kr)]', kO);
dyGxy{1} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]1', kO);
dyGxy{2} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]2', kO);
dyGxy{3} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]3', kO);


uinc = @(X) exp(1i*k0*(cos(theta)*X(:,1) +sin(theta)*X(:,2)));




Uh = fem(m_Gamma, 'P1');
Vh = fem(m_Gamma, 'P0');

G = @(X, Y) femGreenKernel(X, Y, '[H0(kr)]', k0);
dyGxy{1} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]1', k0);
dyGxy{2} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]2', k0);
dyGxy{3} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]3', k0);


S0   = 1i/4 * integral(pts.vtx(~in, :), Gamma, G, Vh);
S0   = S0 -1/(2*pi) .* regularize(pts.vtx(~in, :), Gamma,'[log(r)]', Vh);



D0  = 1i/4 * integral(pts.vtx(~in, :), Gamma, dyGxy, ntimes(Uh));
D0  = D0 -1/(2*pi) .* regularize(pts.vtx(~in, :), Gamma,'grady[log(r)]', ntimes(Uh));




G = @(X, Y) femGreenKernel(X, Y, '[H0(kr)]', kO);
dyGxy{1} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]1', kO);
dyGxy{2} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]2', kO);
dyGxy{3} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]3', kO);


dxGxy{1} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]1', kO);
dxGxy{2} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]2', kO);
dxGxy{3} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]3', kO);


comm_dxGxy{1} = @(X,Y) (alpha(Y) - alpha(X)) .* femGreenKernel(X, Y, 'gradx[H0(kr)]1', kO);
comm_dxGxy{2} = @(X,Y) (alpha(Y) - alpha(X)) .* femGreenKernel(X, Y, 'gradx[H0(kr)]2', kO);
comm_dxGxy{3} = @(X,Y) (alpha(Y) - alpha(X)) .* femGreenKernel(X, Y, 'gradx[H0(kr)]3', kO);



S1   = 1i/4 * integral(pts.vtx(in, :), Gamma, G, Vh);
S1   = S1 -1/(2*pi) .* regularize(pts.vtx(in, :), Gamma,'[log(r)]', Vh);



D1  = 1i/4 * integral(pts.vtx(in, :), Gamma, dyGxy, ntimes(Uh));
D1  = D1 -1/(2*pi) .* regularize(pts.vtx(in, :), Gamma,'grady[log(r)]', ntimes(Uh));


Uh = fem(m_Omega, 'P1');

if (check > 1e-10) && (check2 > 1e-10)


K = (1i/4) * integral(pts.vtx(in, :), Omega, G, Uh);
Kreg = -(1/(2 * pi)) * regularize(pts.vtx(in, :), Omega, '[log(r)]',Uh);


% Beta = diag(beta(m_Omega.vtx));

dim = size(m_Omega.vtx, 1);
Beta = spdiags(beta(m_Omega.vtx), 0, dim, dim); 


% 
aux = uinc(pts.vtx(~in, :));
% uext = aux + (K + Kreg)* (Beta*u);


K = (1i / 4) * integral(pts.vtx(in, :), Omega, dxGxy, grad(Uh));
commK = (1i / 4) * integral(pts.vtx(in, :), Omega, comm_dxGxy, grad(Uh));


reg = - 1 / (2 * pi) * regularize(pts.vtx(in, :), Omega, 'grady[log(r)]', grad(Uh));

m = size(pts.vtx(in, :), 1);
Alpha = spdiags(alpha(pts.vtx(in, :)),0,m,m);



%%
uex = aux - S0 * lambda + D0 * phi;
uin = S1 * lambda - D1 * phi + (K + Kreg)* (Beta*u) + (Alpha * (K - reg) + commK) * u;


    
elseif (check > 1e-10) && (check2 <= 1e-10)

K = (1i/4) * integral(pts.vtx(in, :), Omega, G, Uh);
Kreg = -(1/(2 * pi)) * regularize(pts.vtx(in, :), Omega, '[log(r)]',Uh);


% Beta = diag(beta(m_Omega.vtx));

dim = size(m_Omega.vtx, 1);
Beta = spdiags(beta(m_Omega.vtx), 0, dim, dim); 


% 
aux = uinc(pts.vtx(~in, :));
% uext = aux + (K + Kreg)* (Beta*u);


%%
uex = aux - S0 * lambda + D0 * phi;
uin = S1 * lambda - D1 * phi + (K + Kreg)* (Beta*u);



else
    
    

% 
aux = uinc(pts.vtx(~in, :));
% uext = aux + (K + Kreg)* (Beta*u);


%%
uex = aux - S0 * lambda + D0 * phi;
uin = S1 * lambda - D1 * phi;

end
    
uext = zeros(size(in, 1), 1);

uext(in) = uin;
uext(~in) = uex;


else

beta = @(X) (kx(X).^2) / aO - kO^2;

alpha = @(X) (ax(X) - aO) / aO;
check = norm(beta(m_Omega.vtx));

check2 = norm(alpha(m_Omega.vtx));

% disp("Compute Exterior and Interior Field");

Gamma = dom(m_Gamma, 3);
Omega = dom(m_Omega, 3);


G = @(X, Y) femGreenKernel(X, Y, '[H0(kr)]', kO);
dyGxy{1} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]1', kO);
dyGxy{2} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]2', kO);
dyGxy{3} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]3', kO);


uinc = @(X) exp(1i*k0*(cos(theta)*X(:,1) +sin(theta)*X(:,2)));




Uh = fem(m_Gamma, 'P1');
Vh = fem(m_Gamma, 'P0');

dyGxy{1} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]1', k0);
dyGxy{2} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]2', k0);
dyGxy{3} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]3', k0);









G = @(X, Y) femGreenKernel(X, Y, '[H0(kr)]', kO);
dyGxy{1} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]1', kO);
dyGxy{2} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]2', kO);
dyGxy{3} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]3', kO);


dxGxy{1} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]1', kO);
dxGxy{2} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]2', kO);
dxGxy{3} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]3', kO);


comm_dxGxy{1} = @(X,Y) (alpha(Y) - alpha(X)) .* femGreenKernel(X, Y, 'gradx[H0(kr)]1', kO);
comm_dxGxy{2} = @(X,Y) (alpha(Y) - alpha(X)) .* femGreenKernel(X, Y, 'gradx[H0(kr)]2', kO);
comm_dxGxy{3} = @(X,Y) (alpha(Y) - alpha(X)) .* femGreenKernel(X, Y, 'gradx[H0(kr)]3', kO);



S1   = 1i/4 * integral(pts.vtx(in, :), Gamma, G, Vh);
S1   = S1 -1/(2*pi) .* regularize(pts.vtx(in, :), Gamma,'[log(r)]', Vh);



D1  = 1i/4 * integral(pts.vtx(in, :), Gamma, dyGxy, ntimes(Uh));
D1  = D1 -1/(2*pi) .* regularize(pts.vtx(in, :), Gamma,'grady[log(r)]', ntimes(Uh));


Uh = fem(m_Omega, 'P1');

if (check > 1e-10) && (check2 > 1e-10)


K = (1i/4) * integral(pts.vtx(in, :), Omega, G, Uh);
Kreg = -(1/(2 * pi)) * regularize(pts.vtx(in, :), Omega, '[log(r)]',Uh);


% Beta = diag(beta(m_Omega.vtx));

dim = size(m_Omega.vtx, 1);
Beta = spdiags(beta(m_Omega.vtx), 0, dim, dim); 


% 
% uext = aux + (K + Kreg)* (Beta*u);


K = (1i / 4) * integral(pts.vtx(in, :), Omega, dxGxy, grad(Uh));
commK = (1i / 4) * integral(pts.vtx(in, :), Omega, comm_dxGxy, grad(Uh));


reg = - 1 / (2 * pi) * regularize(pts.vtx(in, :), Omega, 'grady[log(r)]', grad(Uh));

m = size(pts.vtx(in, :), 1);
Alpha = spdiags(alpha(pts.vtx(in, :)),0,m,m);



%%
uin = S1 * lambda - D1 * phi + (K + Kreg)* (Beta*u) + (Alpha * (K - reg) + commK) * u;


    
elseif (check > 1e-10) && (check2 <= 1e-10)

K = (1i/4) * integral(pts.vtx(in, :), Omega, G, Uh);
Kreg = -(1/(2 * pi)) * regularize(pts.vtx(in, :), Omega, '[log(r)]',Uh);


% Beta = diag(beta(m_Omega.vtx));

dim = size(m_Omega.vtx, 1);
Beta = spdiags(beta(m_Omega.vtx), 0, dim, dim); 


% 
% uext = aux + (K + Kreg)* (Beta*u);


%%
uin = S1 * lambda - D1 * phi + (K + Kreg)* (Beta*u);



else
    
    

% 
% uext = aux + (K + Kreg)* (Beta*u);


%%
uin = S1 * lambda - D1 * phi;

end
    

uext = uin;

end