function [Av, Ad, An] = build_2VIOb(m_Gamma, m_Omega, k, beta,aO, ax, grad_ax, alpha, tol)



Omega = dom(m_Omega, 3);
Gamma = dom(m_Gamma, 3);
Vh = fem(m_Omega, 'P1');
gradVh = grad(Vh);

G = @(X, Y) femGreenKernel(X, Y, '[H0(kr)]', k);


dxGxy{1} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]1', k);
dxGxy{2} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]2', k);
dxGxy{3} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]3', k);

dyGxy{1} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]1', k);
dyGxy{2} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]2', k);
dyGxy{3} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]3', k);


 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Assembly of Volume Integral Operator %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PART 1: N(beta u)
K = (1i/4) * integral(Omega, Omega,  Vh, G, Vh, tol);
Kreg = -(1/(2 * pi)) * regularize(Omega,Omega,Vh,'[log(r)]',Vh);

dim = size(m_Omega.vtx, 1);
Beta = spdiags(beta(m_Omega.vtx), 0, dim, dim); 
Alpha = spdiags(alpha(m_Omega.vtx), 0, dim, dim);


%%% grad alpha = grad ax / a1
DxAlpha = spdiags(grad_ax{1}(m_Omega.vtx) / aO, 0, dim, dim); 
DyAlpha = spdiags(grad_ax{2}(m_Omega.vtx) / aO, 0, dim, dim); 

Av = (K + Kreg) * (Beta - k^2 * Alpha);


% PART 2: div N(u grad alpha)
Mv = femLagrangePn(Vh, Omega);


[~,Pv] = Vh.unk;   

Mv = Mv * Pv;

[Xq, Wq, ~] = Omega.qud; 

Wq = spdiags(Wq, 0, size(Wq, 1), size(Wq, 1));

divK = integral(Xq, Omega, dyGxy, Vh, tol);

divK{1} = (1i/4) * divK{1};
divK{2} = (1i/4) * divK{2};

reg{1} = -(1/(2 * pi)) * regularize(Xq, Omega, 'grady[log(r)]1',Vh);
reg{2} = -(1/(2 * pi)) * regularize(Xq, Omega, 'grady[log(r)]2',Vh);


%%% Minus sign: change of variable div_y to div_x
divA = -(Mv' * Wq * (divK{1} + reg{1}) * DxAlpha + Mv' * Wq * (divK{2} + reg{2}) * DyAlpha);

Av = Av - divA;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Assembly of Dirichlet Trace  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Vn = fem(m_Gamma, 'P1');
Ud = fem(m_Gamma, 'P1');

Mv = femLagrangePn(Vn, Gamma); 

[~,Pv] = Vn.unk;   

Mv = Mv * Pv;  

[Xq, Wq, elt2qud] = Gamma.qud;

Wq = spdiags(Wq, 0, size(Wq, 1), size(Wq, 1));


% PART 1: N(beta u)

K = (1i/4) * integral(Xq, Omega, G, Vh, tol);
Kreg = -(1/(2 * pi)) * regularize(Xq, Omega, '[log(r)]',Vh);


dim = size(m_Omega.vtx, 1);
Beta = spdiags(beta(m_Omega.vtx), 0, dim, dim); 
Alpha = spdiags(alpha(m_Omega.vtx), 0, dim, dim);


%%% grad alpha = grad ax / a1
DxAlpha = spdiags(grad_ax{1}(m_Omega.vtx) / aO, 0, dim, dim); 
DyAlpha = spdiags(grad_ax{2}(m_Omega.vtx) / aO, 0, dim, dim); 

Ad = Mv' * Wq * (K + Kreg) * (Beta - k^2 * Alpha);

% PART 2: div N(u grad alpha)
  
divKdir = integral(Xq, Omega, dyGxy, Vh, tol);

divKdir{1} = (1i/4) * divKdir{1};
divKdir{2} = (1i/4) * divKdir{2};

reg{1} = -(1/(2 * pi)) * regularize(Xq, Omega, 'grady[log(r)]1',Vh);
reg{2} = -(1/(2 * pi)) * regularize(Xq, Omega, 'grady[log(r)]2',Vh);


%%% Minus sign: change of variable div_y to div_x
divKdir = -(Mv' * Wq * (divKdir{1} + reg{1}) * DxAlpha + Mv' * Wq * (divKdir{2} + reg{2}) * DyAlpha);

Ad = Ad - divKdir;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Assembly of Neumann Trace  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ud = fem(m_Gamma, 'P1');
nUd = ntimes(Ud);

Mnu = femLagrangePn(nUd, Gamma);
[~,Pnu] = nUd.unk;   

Mnu{1} = Mnu{1} * Pnu; 
Mnu{2} = Mnu{2} * Pnu; 
Mnu{3} = Mnu{3} * Pnu; 


% PART 1: N(beta u)
K = integral(Xq, Omega, dxGxy, Vh, tol);

K{1} = (1i/4) * K{1};
K{2} = (1i/4) * K{2};

clear reg;
reg{1} = -(1/(2 * pi)) * regularize(Xq, Omega, 'grady[log(r)]1',Vh);
reg{2} = -(1/(2 * pi)) * regularize(Xq, Omega, 'grady[log(r)]2',Vh);

An = Mnu{1}' * Wq * (K{1} - reg{1}) + Mnu{2}' * Wq * (K{2} - reg{2});
An = An * (Beta - k^2 * Alpha);


normals = Gamma.qudNrm;

m    = size(normals,1);

Nrm{1} = spdiags(normals(:,1),0,m,m);
Nrm{2} = spdiags(normals(:,2),0,m,m);
Nrm{3} = spdiags(normals(:,3),0,m,m);

% 
% 
% grad_Ax{1} = spdiags(grad_ax{1}(Xq) ./ ax(Xq), 0, size(Wq, 1), size(Wq, 1));
% grad_Ax{2} = spdiags(grad_ax{2}(Xq) ./ ax(Xq), 0, size(Wq, 1), size(Wq, 1));
% grad_Ax{3} = spdiags(grad_ax{3}(Xq) ./ ax(Xq), 0, size(Wq, 1), size(Wq, 1));
% 
% 
% %%% n dot (grad a / a)
% Mult = Nrm{1} * grad_Ax{1} + Nrm{2} * grad_Ax{2} + Nrm{3} * grad_Ax{3};
% 
% Mu = femLagrangePn(Ud, Gamma);
% [~,Pu] = Ud.unk;
% 
% Mu = Mu * Pu;
% 
% 
% K = (1i/4) * integral(Xq, Omega, G, Vh, tol);
% Kreg = -(1/(2 * pi)) * regularize(Xq, Omega, '[log(r)]',Vh);
% 
% %%% Minus sign: - n dot grad a / a
% An1 = -Mu' * Mult * Wq * (K + Kreg) * (Beta - k^2 * Alpha);


% PART 2: div N(u grad alpha)

K1 = (1i / 4) * integral(Xq, Omega, dxGxy{1}, Vh);

K2 = (1i / 4) * integral(Xq, Omega, dxGxy{2}, Vh);

reg{1} = -(1 /(2 * pi)) * regularize(Xq, Omega, 'grady[log(r)]1', Vh);
reg{2} = -(1 /(2 * pi)) * regularize(Xq, Omega, 'grady[log(r)]2', Vh);


Mv = femLagrangePn(grad(Ud), Gamma); 

[~,Pv] = Ud.unk;   

Mv{1} =   Nrm{2} * Mv{1} * Pv; 
Mv{2} =  -Nrm{1} * Mv{2} * Pv; 
Mv{3} = 0*Nrm{3} * Mv{3} * Pv; 
    
%%% n cross (grad G cross grad alpha) u
% I1 = Mv{1}' * Wq * (K1 - reg{1}) * DyAlpha + ... 
%      Mv{2}' * Wq * (K2 - reg{2}) * DxAlpha;

I1 = (Mv{1}' + Mv{2}') * Wq * ((K1 - reg{1}) * DyAlpha - (K2 - reg{2}) * DxAlpha);

 
%%%%%%%%%%%%%%%%%%%%%%%%%%
clear reg;

K = (1i / 4) * integral(Xq, Omega, G, Vh); 

reg = 1 / (2 * pi) * regularize(Xq, Omega, '[log(r)]', Vh);

I2 = Mnu{1}' * Wq * (K - reg) * DxAlpha + ...
Mnu{2}' * Wq * (K - reg) * DyAlpha;

    %%%%%%%%%%%%%%%%%%%%%%% <alpha n dot grad u, phi>_Gamma %%%%%%%%%%%%%%%%%%%
% I3 = integral(Gamma, phih, alpha, psih);

An = An - (I1 - k^2 * I2);    
end