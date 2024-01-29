function [Av, traceA] = build_VIO(m_Gamma, m_Omega, k0, k, beta, tol)


check = norm(beta(m_Omega.vtx));

if check > 1e-10

Omega = dom(m_Omega, 3);
Gamma = dom(m_Gamma, 3);
Vh = fem(m_Omega, 'P1');


G = @(X, Y) femGreenKernel(X, Y, '[H0(kr)]', k);
Gbeta = @(X, Y) G(X, Y) .* beta(X);

Gdiff = @(X, Y) G(X, Y) .* (beta(Y) - beta(X));
% disp('--------------------------------')

% disp('Assembly of Volume Integral Operator');

% tic
K = (1i/4) * integral(Omega, Omega,  Vh, G, Vh, tol);


% Kdiff = (1i/4) * integral(Omega, Omega,  Vh, Gdiff, Vh, tol);
% K = (1i/4) * integral(Omega, Omega,  Vh, Gbeta, Vh, tol);
Kreg = -(1/(2 * pi)) * regularize(Omega,Omega,Vh,'[log(r)]',Vh);

% toc



% if strcmp(typ, 'P0')
%     
%     BetaP0 = integral(Omega, Vh, beta);
%     Mass = integral(Omega, Vh, Vh);
%     
%     Beta = Mass \ BetaP0;
%     
%     Beta = diag(Beta);
%     
% elseif strcmp(typ, 'P1')
% Beta = diag(beta(m_Omega.vtx)); 


dim = size(m_Omega.vtx, 1);
Beta = spdiags(beta(m_Omega.vtx), 0, dim, dim); 
%     
% end
% 
% if strcmp(typ, 'P0')
% A = (K + Kreg) * Beta;
% elseif strcmp(typ, 'P1')


Betah = hmx(double(m_Omega.vtx), double(m_Omega.vtx), Beta, tol);
Av = (K + Kreg) * Betah;
% Av = Kdiff + (K + Beta * Kreg);

% end
% 
% m_Gamma = bnd(m_Omega);

Vn = fem(m_Gamma, 'P0');
Ud = fem(m_Gamma, 'P1');
nUd = ntimes(Ud);

Mv = femLagrangePn(Vn, Gamma); 
Mu = femLagrangePn(Ud, Gamma); 
Mnu = femLagrangePn(nUd, Gamma); 

[~,Pv] = Vn.unk;   
[~,Pu] = Ud.unk;   
[~,Pnu] = nUd.unk;   

Mv = Mv * Pv;  
Mu = Mu * Pu;  
Mnu{1} = Mnu{1} * Pnu; 
Mnu{2} = Mnu{2} * Pnu; 
Mnu{3} = Mnu{3} * Pnu; 


[Xq, Wq, elt2qud] = Gamma.qud; % 81 x 3

Wq = spdiags(Wq, 0, size(Wq, 1), size(Wq, 1));

% disp('--------------------------------')
% disp('Assembly of Dirichlet Trace');
% tic
K = (1i/4) * integral(Gamma, Omega, Vn, G, Vh, tol);
Kreg = -(1/(2 * pi)) * regularize(Gamma, Omega, Vn, '[log(r)]',Vh);


% KKdiff = (1i/4) * integral(Gamma, Omega, Vn, Gdiff, Vh, tol);
% KK = (1i/4) * integral(Gamma, Omega, Vn, Gbeta, Vh, tol);
% KKreg = -(1/(2 * pi)) * regularize(Gamma, Omega, Vn, '[log(r)]',Vh);
% Beta = diag(beta(m_Omega.vtx)); 

% dim = size(m_Omega.vtx, 1);
% Beta = spdiags(beta(m_Omega.vtx), 0, dim, dim); 


% Ad = Mv' * Wq * (K + Kreg) * Beta;

Ad = (K + Kreg) * Betah;

% Ad = KKdiff + (KK + Beta * KKreg);
% toc


% disp('--------------------------------')
% disp('Assembly of Neumann Trace');


dxGxy{1} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]1', k);
dxGxy{2} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]2', k);
dxGxy{3} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]3', k);

dxGxyBeta{1} = @(X,Y) dxGxy{1}(X, Y) .* beta(X);
dxGxyBeta{2} = @(X,Y) dxGxy{2}(X, Y) .* beta(X);
dxGxyBeta{3} = @(X,Y) dxGxy{3}(X, Y) .* beta(X);

dxGxydiff{1} = @(X,Y) dxGxy{1}(X, Y) .* (beta(Y)-beta(X));
dxGxydiff{2} = @(X,Y) dxGxy{2}(X, Y) .* (beta(Y)-beta(X));
dxGxydiff{3} = @(X,Y) dxGxy{3}(X, Y) .* (beta(Y)-beta(X));
% tic

% K = integral(Xq, Omega, dxGxy, Vh, tol);
K = (1i / 4) * integral(Gamma, Omega, ntimes(Ud), dxGxy, Vh, tol);

% K = integral(Xq, Omega, dxGxyBeta, Vh, tol);


% KK = 1i/4 * integral(Gamma, Omega, ntimes(Ud), dxGxyBeta, Vh, tol);
% KKdiff = 1i/4 * integral(Gamma, Omega, ntimes(Ud), dxGxydiff, Vh, tol);

% K{1} = (1i/4) * K{1};
% K{2} = (1i/4) * K{2};

% reg{1} = -(1/(2 * pi)) * regularize(Xq, Omega, 'grady[log(r)]1',Vh);
% reg{2} = -(1/(2 * pi)) * regularize(Xq, Omega, 'grady[log(r)]2',Vh);


% regg = -(1/(2 * pi)) * regularize(Gamma, Omega, ntimes(Ud), 'grady[log(r)]',Vh);

regg = -(1/(2 * pi)) * regularize(Gamma, Omega, ntimes(Ud), 'grady[log(r)]',Vh);



% An = Mnu{1}' * Wq * (K{1} + reg{1}) + Mnu{2}' * Wq * (K{2} + reg{2});
% An = An * Beta;

% An = Mnu{1}' * Wq * (K{1} + reg{1} * Beta) + ...
%      Mnu{2}' * Wq * (K{2} + reg{2} * Beta);

An = (K + regg) * Betah;
% An = KKdiff + (KK + Beta * regg);
 
traceA = [An; Ad];
% toc


else
    
Ng = length(m_Gamma.vtx);
No = length(m_Omega.vtx);
% 
% Av = zeros(No, No);
% Ad = zeros(Ng, No);
% An = zeros(Ng, No);


Mh = hmx();
Mh.typ = 1;
Mh.row = cell(1,4);
Mh.col = cell(1,4);
Mh.chd = cell(1,4);
Mh.dat = {zeros(No,0),zeros(0,No)};

Av = Mh;

    

Mh = hmx();
Mh.typ = 1;
Mh.row = cell(1,4);
Mh.col = cell(1,4);
Mh.chd = cell(1,4);
Mh.dat = {zeros(2*Ng,0),zeros(0,No)};

traceA = Mh;

% disp('--------------------------------')
end