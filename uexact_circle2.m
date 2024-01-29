function [u, uinc] = uexact_circle2(N, mesh, pts, k0, k, theta)


rad = 1;
m_Omega = mesh;
m_Gamma = bnd(mesh);

Omega = dom(mesh, 3); 
Gamma = dom(m_Gamma, 3);

X = m_Gamma.vtx;

cosTheta = (cos(theta) * X(:, 1) + sin(theta) * X(:, 2)) ./ sqrt(X(:, 1).^2 + X(:, 2).^2);

Theta = acos(cosTheta);

Zeta = X(:, 1) + 1i * X(:, 2);

Zeta = Zeta / rad;

Uinc = @(X) exp(1i*k0* (cos(theta) * X(:, 1) + sin(theta) * X(:, 2) ) );



clc
 
u = zeros(size(pts.vtx, 1), 1);

phi = zeros(size(m_Gamma.vtx, 1), 1);

psi = zeros(size(m_Gamma.vtx, 1), 1);

% 
% phi_inc = zeros(size(m_Gamma.vtx, 1), 1);
% 
% psi_inc = zeros(size(m_Gamma.vtx, 1), 1);

for n = -N:N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H0 = besselh(abs(n), 1, k0*rad);
J0 = besselj(abs(n), k0*rad);

Hp0 = besselh(abs(n)+1, 1, k0*rad);
Jp0 = besselj(abs(n)+1, k0*rad);  

Hm0 = besselh(abs(n)-1, 1, k0*rad);
Jm0 = besselj(abs(n)-1, k0*rad);  

dH0 = abs(n) * H0 / (k0*rad) - Hp0;
dJ0 = abs(n) * J0 / (k0*rad) - Jp0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H1 = besselh(abs(n), 1, k*rad);
J1 = besselj(abs(n), k*rad);

Hp1 = besselh(abs(n)+1, 1, k*rad);
Jp1 = besselj(abs(n)+1, k*rad);  

Hm1 = besselh(abs(n)-1, 1, k*rad);
Jm1 = besselj(abs(n)-1, k*rad);  

dH1 = abs(n) * H1 / (k*rad) - Hp1;
dJ1 = abs(n) * J1 / (k*rad) - Jp1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beta = (-1)^(abs(n) +1) * (k0 * dJ0 * H0 - k * J0 * dH0) / ... 
                          (k0 * dH0 * J1 - k * H0 * dJ1);
                      
phi = phi + beta * J1 * cos(n * Theta);

psi = psi + beta * k * dJ1 * cos(n * Theta);


end

% phi = real(phi);
% psi = real(psi);


% phi_inc = real(phi_inc);
% psi_inc = real(psi_inc);
% clc
% 
% ref = Uinc(X);
% ref = ref(2:end);
% error = norm(uinc(2:end) - ref);

% disp(error)

Gamma = dom(m_Gamma, 3);
Vh = fem(m_Gamma, 'P0');
Uh = fem(m_Gamma, 'P1');

G = @(X, Y) femGreenKernel(X, Y, '[H0(kr)]', k);
dyGxy{1} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]1', k);
dyGxy{2} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]2', k);
dyGxy{3} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]3', k);


S1   = 1i/4 * integral(pts.vtx, Gamma, G, Uh);
S1   = S1 -1/(2*pi) .* regularize(pts.vtx, Gamma,'[log(r)]', Uh);



D1  = 1i/4 * integral(pts.vtx, Gamma, dyGxy, ntimes(Uh));
D1  = D1 -1/(2*pi) .* regularize(pts.vtx, Gamma,'grady[log(r)]', ntimes(Uh));

u = S1 * psi - D1 * phi;
% uinc = S1 * psi_inc - D1 * phi_inc;
uinc = Uinc(pts.vtx);

% Wh = fem(m_Omega, 'P1');
% [~,I1,I2] = intersect(Uh.unk,Wh.unk,'rows','stable');
% 
% u(I2) = phi;

% disp("Ok");

end