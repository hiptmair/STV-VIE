addpath('openMsh');
addpath('openDom');
addpath('openFem');
addpath('openHmx');
addpath('miscellaneous');  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% u - div N (alpha grad u) - N (beta u) %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
 
N = 70;
k0 = 1;
k = 2;
omega = 4;
theta = -pi/4;
a0 = 1;

for type = 3

for Const = [1.5]

rad = 1;


f = @(X) 0 * X(:, 1);

kx = @(X) k * ones(size(X(:, 1)));

% kx = @(X) k + 2 * eta3(rad, X);
% kx = @(X) k + 2*eta(rad*sqrt((X(:, 1)).^2 + (X(:, 2)).^2));
% 
% ax = @(X) Const + eta2(rad, X);
% 
% grad_ax{1} = @(X) etap2(rad, X, 1);
% 
% grad_ax{2} = @(X) etap2(rad, X, 2);
%               
% grad_ax{3} = @(X) etap2(rad, X, 3);


% 
% ax = @(X) Const + 1 - X(:, 1).^2 - X(:, 2).^2;
% 
% grad_ax{1} = @(X) -2 * X(:, 1);
% grad_ax{2} = @(X) -2 * X(:, 2);
% grad_ax{3} = @(X) -2 * X(:, 3);



ax = @(X) Const * ones(size(X(:, 1)));

grad_ax{1} = @(X) zeros(size(X(:, 1)));

grad_ax{2} = @(X) zeros(size(X(:, 1)));
              
grad_ax{3} = @(X) zeros(size(X(:, 1)));

% % 
% rad2 = 1;
% % 
% 
% ax = @(X) Const + 0.5 * (a0 - Const) * eta3(rad2, X);
% 
% grad_ax{1} = @(X) 0.5 * (a0 - Const) * etap3(rad2, X, 1);
% 
% grad_ax{2} = @(X) 0.5 * (a0 - Const) * etap3(rad2, X, 2);
%               
% grad_ax{3} = @(X) 0.5 * (a0 - Const) * etap3(rad2, X, 3);


pts = mshSquare(100, [4 4]);
% % 
% ax = @(X) Const + 0.5 - 0.25 * (X(:, 1).^2 + X(:, 2).^2);
% 
% grad_ax{1} = @(X) -0.25 * 2 * X(:, 1);
% grad_ax{2} = @(X) -0.25 * 2 * X(:, 2);
% grad_ax{3} = @(X) 0 * X(:, 3);
% 
% 
% ax = @(X) Const + 0.5 * cos(X(:, 1) * pi);
% 
% grad_ax{1} = @(X) - 0.5 * pi * sin(X(:, 1) * pi);
% grad_ax{2} = @(X) 0 * X(:, 2);
% grad_ax{3} = @(X) 0 * X(:, 3);

% Nref = 400;
% Ns = [100 200];
% count = 1;
typ = 'P1';

% m = mshDisk(Nref, 1);


% mbnd = bnd(m);
% m_Gamma = mbnd;
tol = 1e-3;
% in = inpolygon(pts.vtx(:, 1), pts.vtx(:, 2), m_Gamma.vtx(:, 1), m_Gamma.vtx(:, 2));


% m = mshDisk(N, 1);  
% m = mshSquare(N, [2 2]);

m = mshSquareSector(0, type);

meshes{1} = m;
for j = 1:2
%    meshes{j+1} = meshes{j}.refine; 
   I = (1:size(meshes{j}.elt, 1))';
   meshes{j+1} = mshMidpoint(meshes{j}, I); 
end

mref = meshes{end};

mbnd_ref = bnd(mref);

if type ==1
xS = [-1 -1;1 -1;0 0;1 0;1 1;-1 1;-1 -1];

elseif type == 2
xS = [-1 -1;0 -1;0 0;1 0;1 1;-1 1;-1 -1];

elseif type == 3
xS = [-1 -1;1 -1;1 1;-1 1;-1 -1];

end
in = inpolygon(m.vtx(:, 1), m.vtx(:, 2), xS(:, 1), xS(:, 2));

tic
% [uh_ref, ~] = solve_VIE2(mref, m, k0, k, kx, ax, typ, tol, theta);

[uh_ref, ~] = solve_STF_VIEb(mbnd_ref, mref, m, in, k0, k, 0, kx, ax, grad_ax, tol, theta);

toc

errors = zeros(5, length(meshes)-1);
errorsH1 = zeros(5, length(meshes)-1);
steps = zeros(5, length(meshes)-1);


for count = 1:(length(meshes)-1)

m = meshes{count};
Omega = dom(m, 3);
Uh = fem(m, 'P1');
   

tic
[uh, ~] = solve_VIE2(m, pts, k0, k, kx, ax, typ, tol, theta);

errors(1,count) = compute_errors3(m, mref, uh, uh_ref, 'L2');
errorsH1(1,count) = compute_errors3(m, mref, uh, uh_ref, 'H1');
toc

mbnd = bnd(m);

% in = inpolygon(pts.vtx(:, 1), pts.vtx(:, 2), mbnd.vtx(:, 1), mbnd.vtx(:, 2));

in = inpolygon(pts.vtx(:, 1), pts.vtx(:, 2), xS(:, 1), xS(:, 2));

tic
[uhh, phih, ~] = solve_CVIEb(mbnd, m, pts, k0, k, kx, 1.0, ax, grad_ax, typ, tol, theta);


errors(2,count) = compute_errors3(m, mref, uhh, uh_ref, 'L2');
errorsH1(2,count) = compute_errors3(m, mref, uhh, uh_ref, 'H1');
toc

tic
[u2h, ~] = solve_STF_VIEb(mbnd, m, pts, in, k0, k, 0, kx, ax, grad_ax, tol, theta);

errors(3,count) = compute_errors3(m, mref, u2h, uh_ref, 'L2');
errorsH1(3,count) = compute_errors3(m, mref, u2h, uh_ref, 'H1');
toc


tic
[u3h, ~] = solve_2STF_VIEb(mbnd, m, pts, in, k0, k, 0, kx, ax, grad_ax, tol, theta);

errors(4,count) = compute_errors3(m, mref, u3h, uh_ref, 'L2');
errorsH1(4,count) = compute_errors3(m, mref, u3h, uh_ref, 'H1');
toc

tic
uFEM = solve_FEM_BEM(m, k0, kx, theta, k, ax, f);

errors(5,count) = compute_errors3(m, mref, uFEM, uh_ref, 'L2');
errorsH1(5,count) = compute_errors3(m, mref, uFEM, uh_ref, 'H1');
toc


step = m.ndv;
steps(1,count) = sqrt(max(step));

disp('-----------------------------------------')
end

% [order, orderH1] = convergence_rates(errors, errorsH1, meshes, steps);
% 
% if type == 1
% save(['resultsC/Zshaped_', num2str(Const), '.mat'], 'errors', 'order', 'errorsH1', 'orderH1', 'steps')
% 
% elseif type == 2
% save(['resultsC/Lshaped_', num2str(Const), '.mat'], 'errors', 'order', 'errorsH1', 'orderH1', 'steps')
% 
% elseif type == 3
% save(['resultsC/Square_', num2str(Const), '.mat'], 'errors', 'order', 'errorsH1', 'orderH1', 'steps')
% end

disp('-----------------------------------------')
disp('-----------------------------------------')
end

end

%%%%%%%%%%%%%%%%
%%





% convergence_plot(steps, errors, errorsH1);



%%




% 
% Vh = fem(mref, 'P1');
% 
% figure
% graph(Vh, real(uh_ref));hold on;
% title('STF-VIE','Interpreter', 'latex', 'FontSize', 18);
% xlim([-1 1]);
% ylim([-1 1])
% axis equal
% colormap 'jet'
% colorbar 
% caxis([-1.2 1.2]);
