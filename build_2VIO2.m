function [Av, Ad, An] = build_2VIO2(m_Gamma, m_Omega, k, beta, a, alpha, tol)


check = norm(beta(m_Omega.vtx));

check2 = norm(alpha(m_Omega.vtx));

if (check > 1e-10) && (check2 > 1e-10)

Omega = dom(m_Omega, 3);
Gamma = dom(m_Gamma, 3);
Vh = fem(m_Omega, 'P1');
gradVh = grad(Vh);
Gu = gradVh.uqm(Omega);


G = @(X, Y) femGreenKernel(X, Y, '[H0(kr)]', k);
comm_G = @(X, Y) (alpha(Y) - alpha(X)) .* femGreenKernel(X, Y, '[H0(kr)]', k);


dxGxy{1} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]1', k);
dxGxy{2} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]2', k);
dxGxy{3} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]3', k);


comm_dxGxy{1} = @(X,Y) (alpha(Y) - alpha(X)) .* femGreenKernel(X, Y, 'gradx[H0(kr)]1', k);
comm_dxGxy{2} = @(X,Y) (alpha(Y) - alpha(X)) .* femGreenKernel(X, Y, 'gradx[H0(kr)]2', k);
comm_dxGxy{3} = @(X,Y) (alpha(Y) - alpha(X)) .* femGreenKernel(X, Y, 'gradx[H0(kr)]3', k);
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Assembly of Volume Integral Operator %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PART 1: N(beta u)
K = (1i/4) * integral(Omega, Omega,  Vh, G, Vh, tol);
Kreg = -(1/(2 * pi)) * regularize(Omega,Omega,Vh,'[log(r)]',Vh);






dim = size(m_Omega.vtx, 1);
Beta = spdiags(beta(m_Omega.vtx), 0, dim, dim); 

Av = (K + Kreg) * Beta;


% PART 2: div N(alpha grad u)

[Xq, Wq, elt2qud] = Omega.qud;
Wq = spdiags(Wq, 0, size(Wq, 1), size(Wq, 1));


K = (1i / 4) * integral(Xq, Omega, dxGxy, grad(Vh), tol);
commK = (1i / 4) * integral(Xq, Omega, comm_dxGxy, grad(Vh), tol);


reg = - 1 / (2 * pi) * regularize(Xq, Omega, 'grady[log(r)]', grad(Vh));

m = size(Xq, 1);
Alpha = spdiags(alpha(Xq),0,m,m);

Mv = femLagrangePn(Vh, Omega); 

[~,Pv] = Vh.unk;  

Mv = Mv * Pv;


Av = Av + Mv' * Wq * (Alpha * (K - reg) + commK);

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

Ad = Mv' * Wq * (K + Kreg) * Beta;

% PART 2: div N(alpha grad u)
K = (1i / 4) * integral(Xq, Omega, dxGxy, grad(Vh), tol);

commK = (1i / 4) * integral(Xq, Omega, comm_dxGxy, grad(Vh), tol);


reg = - 1 / (2 * pi) * regularize(Xq, Omega, 'grady[log(r)]', grad(Vh));

m = size(Xq, 1);
Alpha = spdiags(alpha(Xq),0,m,m);

if sum(sum(abs(Alpha))) < 1e-12
    
    Ad = Ad + Mv' * Wq * commK;
else
    
    Ad = Ad + Mv' * Wq * (Alpha * (K - reg) + commK);
end

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

An = Mnu{1}' * Wq * (K{1} + reg{1}) + Mnu{2}' * Wq * (K{2} + reg{2});
An = An * Beta;



% PART 2: div N(alpha grad u)

% Mu = femLagrangePn(Ud, Gamma);
% [~,Pu] = Ud.unk;
% 
% Mu = Mu * Pu;



% int_Gamma grad(phi) dot alpha(x)tau(x) int_Omega gradG x gradU
%+int_Gamma grad(phi) dot tau(x) int_Omega comm_gradG x gradU

if sum(sum(abs(Alpha))) < 1e-12
    

    K1 = (1i / 4) * integral(Xq, Omega, comm_dxGxy{1}, grad(Vh, 2));

    K2 = (1i / 4) * integral(Xq, Omega, comm_dxGxy{2}, grad(Vh, 1));

    tau = Gamma.qudTgt;
    normals = Gamma.qudNrm;

    m    = size(tau,1);
    T{1} = spdiags(tau(:,1),0,m,m);
    T{2} = spdiags(tau(:,2),0,m,m);


    Nrm{1} = spdiags(normals(:,1),0,m,m);
    Nrm{2} = spdiags(normals(:,2),0,m,m);
    Nrm{3} = spdiags(normals(:,3),0,m,m);

    Mv = femLagrangePn(grad(Ud), Gamma); 

    [~,Pv] = Ud.unk;   

    Mv{1} =   Nrm{2} * Mv{1} * Pv; 
    Mv{2} =  -Nrm{1} * Mv{2} * Pv; 
    Mv{3} = 0*Nrm{3} * Mv{3} * Pv; 
    
    
    I1 = -(Mv{1}' + Mv{2}') * Wq * (K1 - K2);
else
    
    K1 = (1i / 4) * integral(Xq, Omega, dxGxy{1}, grad(Vh, 2));

    K2 = (1i / 4) * integral(Xq, Omega, dxGxy{2}, grad(Vh, 1));

    reg = -(1 /(2 * pi)) * regularize(Xq, Omega, 'grady[log(r)]', xgrad(Vh));

    K = K1 - K2 - reg;
    


    K1 = (1i / 4) * integral(Xq, Omega, comm_dxGxy{1}, grad(Vh, 2));

    K2 = (1i / 4) * integral(Xq, Omega, comm_dxGxy{2}, grad(Vh, 1));

    tau = Gamma.qudTgt;
    normals = Gamma.qudNrm;

    m    = size(tau,1);
    T{1} = spdiags(tau(:,1),0,m,m);
    T{2} = spdiags(tau(:,2),0,m,m);

    Nrm{1} = spdiags(normals(:,1),0,m,m);
    Nrm{2} = spdiags(normals(:,2),0,m,m);
    Nrm{3} = spdiags(normals(:,3),0,m,m);

    Mv = femLagrangePn(grad(Ud), Gamma); 

    [~,Pv] = Ud.unk;   

    Mv{1} =   Nrm{2} * Mv{1} * Pv; 
    Mv{2} =  -Nrm{1} * Mv{2} * Pv; 
    Mv{3} = 0*Nrm{3} * Mv{3} * Pv; 
    
    
    I1 = -(Mv{1}' + Mv{2}') * Wq * (Alpha * K + K1 - K2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
K = integral(Xq, Omega, G, grad(Vh)); 

if sum(sum(abs(Alpha))) < 1e-12
    
    K = (1i / 4) * integral(Gamma, Omega, ntimes(Ud), comm_G, grad(Vh));

    I2 = K;
else
    
    K{1} = (1i / 4) * K{1};
    K{2} = (1i / 4) * K{2};

    reg1 = 1 / (2 * pi) * regularize(Xq, Omega, '[log(r)]', grad(Vh, 1));
    reg2 = 1 / (2 * pi) * regularize(Xq, Omega, '[log(r)]', grad(Vh, 2));
    
    I2 = Mnu{1}' * Wq * (Alpha * (K{1} - reg1)) + ...
    Mnu{2}' * Wq * (Alpha * (K{2} - reg2));

    K = (1i / 4) * integral(Gamma, Omega, ntimes(Ud), comm_G, grad(Vh));

    I2 = I2 + K;
end
%%%%%%%%%%%%%%%%%%%%%%% <alpha n dot grad u, phi>_Gamma %%%%%%%%%%%%%%%%%%%
% I3 = integral(Gamma, phih, alpha, psih);

An = An - I1 - k^2 * I2;


elseif (check > 1e-10) && (check2 <= 1e-10)
    

Omega = dom(m_Omega, 3);
Gamma = dom(m_Gamma, 3);
Vh = fem(m_Omega, 'P1');
gradVh = grad(Vh);

G = @(X, Y) femGreenKernel(X, Y, '[H0(kr)]', k);


dxGxy{1} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]1', k);
dxGxy{2} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]2', k);
dxGxy{3} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]3', k);
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Assembly of Volume Integral Operator %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PART 1: N(beta u)
K = (1i/4) * integral(Omega, Omega,  Vh, G, Vh, tol);
Kreg = -(1/(2 * pi)) * regularize(Omega,Omega,Vh,'[log(r)]',Vh);






dim = size(m_Omega.vtx, 1);
Beta = spdiags(beta(m_Omega.vtx), 0, dim, dim); 

Av = (K + Kreg) * Beta;
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

Ad = Mv' * Wq * (K + Kreg) * Beta;


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

An = Mnu{1}' * Wq * (K{1} + reg{1}) + Mnu{2}' * Wq * (K{2} + reg{2});
An = An * Beta;


elseif (check <= 1e-10) && (check2 > 1e-10)

Omega = dom(m_Omega, 3);
Gamma = dom(m_Gamma, 3);
Vh = fem(m_Omega, 'P1');
gradVh = grad(Vh);
Gu = gradVh.uqm(Omega);


G = @(X, Y) femGreenKernel(X, Y, '[H0(kr)]', k);
comm_G = @(X, Y) (alpha(Y) - alpha(X)) .* femGreenKernel(X, Y, '[H0(kr)]', k);


dxGxy{1} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]1', k);
dxGxy{2} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]2', k);
dxGxy{3} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]3', k);


comm_dxGxy{1} = @(X,Y) (alpha(Y) - alpha(X)) .* femGreenKernel(X, Y, 'gradx[H0(kr)]1', k);
comm_dxGxy{2} = @(X,Y) (alpha(Y) - alpha(X)) .* femGreenKernel(X, Y, 'gradx[H0(kr)]2', k);
comm_dxGxy{3} = @(X,Y) (alpha(Y) - alpha(X)) .* femGreenKernel(X, Y, 'gradx[H0(kr)]3', k);
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Assembly of Volume Integral Operator %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PART 2: div N(alpha grad u)

[Xq, Wq, elt2qud] = Omega.qud;
Wq = spdiags(Wq, 0, size(Wq, 1), size(Wq, 1));


K = (1i / 4) * integral(Xq, Omega, dxGxy, grad(Vh), tol);
commK = (1i / 4) * integral(Xq, Omega, comm_dxGxy, grad(Vh), tol);


reg = - 1 / (2 * pi) * regularize(Xq, Omega, 'grady[log(r)]', grad(Vh));

m = size(Xq, 1);
Alpha = spdiags(alpha(Xq),0,m,m);

Mv = femLagrangePn(Vh, Omega); 

[~,Pv] = Vh.unk;  

Mv = Mv * Pv;


Av = Mv' * Wq * (Alpha * (K - reg) + commK);

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

% PART 2: div N(alpha grad u)
K = (1i / 4) * integral(Xq, Omega, dxGxy, grad(Vh), tol);

commK = (1i / 4) * integral(Xq, Omega, comm_dxGxy, grad(Vh), tol);


reg = - 1 / (2 * pi) * regularize(Xq, Omega, 'grady[log(r)]', grad(Vh));

m = size(Xq, 1);
Alpha = spdiags(alpha(Xq),0,m,m);

Ad = Mv' * Wq * (Alpha * (K - reg) + commK);

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

% PART 2: div N(alpha grad u)

% Mu = femLagrangePn(Ud, Gamma);
% [~,Pu] = Ud.unk;
% 
% Mu = Mu * Pu;

K1 = (1i / 4) * integral(Xq, Omega, dxGxy{1}, grad(Vh, 2));

K2 = (1i / 4) * integral(Xq, Omega, dxGxy{2}, grad(Vh, 1));

reg = -(1 /(2 * pi)) * regularize(Xq, Omega, 'grady[log(r)]', xgrad(Vh));

K = K1 - K2 - reg;

K1 = (1i / 4) * integral(Xq, Omega, comm_dxGxy{1}, grad(Vh, 2));

K2 = (1i / 4) * integral(Xq, Omega, comm_dxGxy{2}, grad(Vh, 1));

tau = Gamma.qudTgt;
normals = Gamma.qudNrm;

m    = size(tau,1);
T{1} = spdiags(tau(:,1),0,m,m);
T{2} = spdiags(tau(:,2),0,m,m);
T{3} = spdiags(tau(:,3),0,m,m);


Nrm{1} = spdiags(normals(:,1),0,m,m);
Nrm{2} = spdiags(normals(:,2),0,m,m);
Nrm{3} = spdiags(normals(:,3),0,m,m);

Mv = femLagrangePn(grad(Ud), Gamma); 

[~,Pv] = Ud.unk;   
  
Mv{1} =   Nrm{2} * Mv{1} * Pv; 
Mv{2} =  -Nrm{1} * Mv{2} * Pv; 
Mv{3} = 0*Nrm{3} * Mv{3} * Pv; 

% int_Gamma grad(phi) dot alpha(x)tau(x) int_Omega gradG x gradU
%+int_Gamma grad(phi) dot tau(x) int_Omega comm_gradG x gradU
I1 = -(Mv{1}' + Mv{2}') * Wq * (Alpha * K + K1 - K2);


K = integral(Xq, Omega, G, grad(Vh));
K{1} = (1i / 4) * K{1};
K{2} = (1i / 4) * K{2};

reg1 = 1 / (2 * pi) * regularize(Xq, Omega, '[log(r)]', grad(Vh, 1));
reg2 = 1 / (2 * pi) * regularize(Xq, Omega, '[log(r)]', grad(Vh, 2));
 

% int_Gamma grad(phi) dot alpha(x)tau(x) int_Omega gradG x gradU
%+int_Gamma grad(phi) dot tau(x) int_Omega comm_gradG x gradU
I2 = Mnu{1}' * Wq * (Alpha * (K{1} - reg1)) + ...
     Mnu{2}' * Wq * (Alpha * (K{2} - reg2));


K = (1i / 4) * integral(Gamma, Omega, ntimes(Ud), comm_G, grad(Vh));

I2 = I2 + K;


%%%%%%%%%%%%%%%%%%%%%%% <alpha n dot grad u, phi>_Gamma %%%%%%%%%%%%%%%%%%%
% I3 = integral(Gamma, phih, alpha, psih);

An = -I1 - k^2 * I2;



elseif (check <= 1e-10) && (check2 <= 1e-10)
    
Ng = length(m_Gamma.vtx);
No = length(m_Omega.vtx);

Av = zeros(No, No);
Ad = zeros(Ng, No);
An = zeros(Ng, No);

    
end