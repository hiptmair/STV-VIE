function [u, uinc] = uexact_circle(N, mesh, pts, k0, k, theta)


rad = 1;
m_Omega = mesh;
m_Gamma = bnd(mesh);

Omega = dom(mesh, 3); 
Gamma = dom(m_Gamma, 3);

X = m_Gamma.vtx;

cosTheta = (cos(theta) * X(:, 1) + sin(theta) * X(:, 2)) ./ sqrt(X(:, 1).^2 + X(:, 2).^2);

Theta = acos(cosTheta);

Uinc = @(X) exp(1i*k0* (cos(theta) * X(:, 1) + sin(theta) * X(:, 2) ) );



clc
 
u = zeros(size(pts.vtx, 1), 1);

phi = zeros(size(m_Gamma.vtx, 1), 1);

psi = zeros(size(m_Gamma.vtx, 1), 1);


phi_inc = zeros(size(m_Gamma.vtx, 1), 1);

psi_inc = zeros(size(m_Gamma.vtx, 1), 1);

for n = -N:N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = besselh(abs(n), 1, k0*rad);
J = besselj(abs(n), k0*rad);

Hp = besselh(abs(n)+1, 1, k0*rad);
Jp = besselj(abs(n)+1, k0*rad);  

Hm = besselh(abs(n)-1, 1, k0*rad);
Jm = besselj(abs(n)-1, k0*rad);  

dH = abs(n) * H / (k0*rad) - Hp;
dJ = abs(n) * J / (k0*rad) - Jp;


V0 = 1i * rad^2 * pi^2 * H * J;

K0 = -1i * k0 * rad^2 * pi^2 * 0.5 * (H * dJ + dH * J);

Kp0 = 1i * k0 * rad^2 * pi^2 * 0.5 * (H * dJ + dH * J);

W0 = -1i * k0^2 * rad^2 * pi^2 * dH * dJ;

A0 = [-K0 V0;-W0 Kp0];
% 
% uinc = (1i)^n * J .* cos(n*Theta);
% 
% duinc = 1i * k0 * cos(Theta) .* uinc;

uinc = 2 * pi * rad * (1i)^n * J;

duinc = 2 * pi * rad * k0 * (1i)^n * dJ;

uinc = [uinc; duinc];


phi_inc = phi_inc + (1i)^n * J .* cos(n*Theta);

psi_inc = psi_inc + k0 * (1i)^n * dJ .* cos(n*Theta);
% 
% phi_inc = phi_inc + (1i)^n * J .* exp(1i*n*Theta);
% 
% psi_inc = psi_inc + 1i * k0 * (1i)^n * dJ .* exp(1i*n*Theta);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = besselh(abs(n), 1, k*rad);
J = besselj(abs(n), k*rad);

Hp = besselh(abs(n)+1, 1, k*rad);
Jp = besselj(abs(n)+1, k*rad);  

Hm = besselh(abs(n)-1, 1, k*rad);
Jm = besselj(abs(n)-1, k*rad);  

dH = abs(n) * H / (k*rad) - Hp;
dJ = abs(n) * J / (k*rad) - Jp;


V1 = 1i * rad^2 * pi^2 * H * J;

K1 = -1i * k * rad^2 * pi^2 * 0.5 * (H * dJ + dH * J);

Kp1 = 1i * k * rad^2 * pi^2 * 0.5 * (H * dJ + dH * J);

W1 = -1i * k^2 * rad^2 * pi^2 * dH * dJ;

A1 = [-K1 V1;-W1 Kp1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sol = (A0 + A1) \ uinc;    
% 
phi = phi + sol(1) .* cos(n*Theta);

psi = psi + sol(2) .* cos(n*Theta);

% 
% phi = phi + sol(1) .* exp(1i*n*Theta);
% 
% psi = psi + sol(2) .* exp(1i*n*Theta);





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