function [A0, A1] = build_STFbb(mesh, k0, k, aO, a, grad_a, tol)


Gamma = dom(mesh, 3);
Vh = fem(mesh, 'P0');
Uh = fem(mesh, 'P1');
nUh = ntimes(Uh);
curlUh = nxgrad(Uh);
%%


Gxy = @(X,Y) femGreenKernel(X, Y, '[H0(kr)]', k0);


dyGxy{1} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]1', k0);
dyGxy{2} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]2', k0);
dyGxy{3} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]3', k0);


dxGxy{1} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]1', k0);
dxGxy{2} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]2', k0);
dxGxy{3} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]3', k0);
%%
% disp('--------------------------------')
% disp('Assembly of Integral Operator A0');
% tic

Mv = femLagrangePn(Vh, Gamma);
Mu = femLagrangePn(Uh, Gamma);
Mnu = femLagrangePn(nUh, Gamma);
Mcu = femLagrangePn(curlUh, Gamma);
[~,Pv] = Vh.unk; 
[~,Pu] = Uh.unk; 
[~,Pnu] = nUh.unk; 
[~,Pcu] = curlUh.unk; 

Mv = Mv * Pv;  
Mu = Mu * Pu;  
Mnu{1} = Mnu{1} * Pnu; 
Mnu{2} = Mnu{2} * Pnu; 
Mnu{3} = Mnu{3} * Pnu; 

Mcu{1} = Mcu{1} * Pcu; 
Mcu{2} = Mcu{2} * Pcu; 
Mcu{3} = Mcu{3} * Pcu; 


[Xq, Wq, elt2qud] = Gamma.qud;

Wq = spdiags(Wq, 0, size(Wq, 1), size(Wq, 1));
Ax = spdiags(a(Xq), 0, size(Wq, 1), size(Wq, 1));


% mAx = spdiags(1./a(Xq), 0, size(Wq, 1), size(Wq, 1));
mAx = spdiags(1./a(mesh.vtx), 0, size(mesh.vtx, 1), size(mesh.vtx, 1));


V = (1i/4) * integral(Xq, Gamma, Gxy, Vh);
V = V - 1/2/pi * regularize(Xq, Gamma,'[log(r)]',Vh);

V = (Mv' * Wq) * (Ax * V); 


K = (1i/4) * integral(Xq, Gamma, dyGxy, ntimes(Uh));
K = K - 1/(2*pi) * regularize(Xq, Gamma, 'grady[log(r)]',ntimes(Uh));

K = Mv' * Wq * Ax * K * mAx; 

% K = (1i/4) * integral(Gamma, Gamma, Vh, dyGxy, ntimes(Uh));
% K = K - 1/(2*pi) * regularize(Gamma, Gamma, Vh,'grady[log(r)]',ntimes(Uh));
% 
% K = K + (1i/4) * integral(Gamma, Gamma, Vh, comm_dyGxy, ntimes(Uh));


Kp = (1i/4) * integral(Gamma, Gamma, ntimes(Uh), dxGxy, Vh);
Kp = Kp - 1/(2*pi) * regularize(Gamma, Gamma, Vh,'grady[log(r)]',ntimes(Uh)).';

W = 1i/4 * (k0^2 * integral(Gamma, Gamma, ntimes(Uh),Gxy,ntimes(Uh)) ...
    - integral(Gamma, Gamma, nxgrad(Uh), Gxy, nxgrad(Uh)));
W = W -1/(2*pi) * (k0^2 * regularize(Gamma, Gamma, ntimes(Uh),'[log(r)]',ntimes(Uh)) ...
    - regularize(Gamma, Gamma, nxgrad(Uh),'[log(r)]',nxgrad(Uh)));

W = W * mAx;

% comm_W = 1i/4 * (k0^2 * integral(Gamma, Gamma, ntimes(Uh),comm_Gxy,ntimes(Uh)) ...
%     - integral(Gamma, Gamma, nxgrad(Uh), comm_Gxy, nxgrad(Uh)));
% 
% 
% W1 = integral(Xq, Gamma, Gxy, ntimes(Uh));
% W2 = integral(Xq, Gamma, Gxy, nxgrad(Uh));
% 
% reg1 = regularize(Xq, Gamma,'[log(r)]',ntimes(Uh));
% reg2 = regularize(Xq, Gamma,'[log(r)]',nxgrad(Uh));
% 
% W1{1} = 1i/4 * k0^2 * W1{1} -1/(2*pi) * k0^2 * reg1{1};
% W1{2} = 1i/4 * k0^2 * W1{2} -1/(2*pi) * k0^2 * reg1{2};
% W1{3} = 1i/4 * k0^2 * W1{3} -1/(2*pi) * k0^2 * reg1{3};
% 
% W2{1} =-1i/4 * W2{1} +1/(2*pi) * reg2{1};
% W2{2} =-1i/4 * W2{2} +1/(2*pi) * reg2{2};
% W2{3} =-1i/4 * W2{3} +1/(2*pi) * reg2{3};
% 
% W1 = Mnu{1}' * Wq * mAx * W1{1} + ...
%      Mnu{2}' * Wq * mAx * W1{2}+ ...
%      Mnu{3}' * Wq * mAx * W1{3};
% 
% 
% W2 = Mcu{1}' * Wq * mAx * W2{1} + ...
%      Mcu{2}' * Wq * mAx * W2{2}+ ...
%      Mcu{3}' * Wq * mAx * W2{3};

% W = W1 + W2 + comm_W;

A0 = [-K V;-W Kp];


%%
Gxy = @(X,Y) femGreenKernel(X, Y, '[H0(kr)]', k);
dyGxy{1} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]1', k);
dyGxy{2} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]2', k);
dyGxy{3} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]3', k);

dxGxy{1} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]1', k);
dxGxy{2} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]2', k);
dxGxy{3} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]3', k);


%%
% disp('--------------------------------')
% disp('Assembly of Integral Operator A1');
% tic

V = (1i/4) * integral(Gamma, Gamma, Vh, Gxy, Vh);
V = V - 1/2/pi * regularize(Gamma, Gamma, Vh,'[log(r)]',Vh);

K = (1i/4) * integral(Gamma, Gamma, Vh, dyGxy, ntimes(Uh));
K = K - 1/(2*pi) * regularize(Gamma, Gamma, Vh,'grady[log(r)]',ntimes(Uh));

Kp = (1i/4) * integral(Gamma, Gamma, ntimes(Uh), dxGxy, Vh);
Kp = Kp - 1/(2*pi) * regularize(Gamma, Gamma, Vh,'grady[log(r)]',ntimes(Uh)).';

W = 1i/4 * (k^2 * integral(Gamma, Gamma, ntimes(Uh),Gxy,ntimes(Uh)) ...
    - integral(Gamma, Gamma, nxgrad(Uh), Gxy, nxgrad(Uh)));
W = W -1/(2*pi) * (k^2 * regularize(Gamma, Gamma, ntimes(Uh),'[log(r)]',ntimes(Uh)) ...
    - regularize(Gamma, Gamma, nxgrad(Uh),'[log(r)]',nxgrad(Uh)));
% toc

% 
% Axm1_gradAx{1} = @(X) grad_a{1}(Xq) ./ a(Xq);
% Axm1_gradAx{2} = @(X) grad_a{2}(Xq) ./ a(Xq);
% Axm1_gradAx{3} = @(X) grad_a{3}(Xq) ./ a(Xq);

% dn_alpha = integral(Gamma, ntimes(Uh), Axm1_gradAx, Uh);

A1 = [-K V;-W Kp];


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%% Operator T_Gamma %%%%%%%%%%%%%
% normals = Gamma.qudNrm;
% 
% m    = size(normals,1);
% 
% Nrm{1} = spdiags(normals(:,1),0,m,m);
% Nrm{2} = spdiags(normals(:,2),0,m,m);
% Nrm{3} = spdiags(normals(:,3),0,m,m);
% 
% grad_Ax{1} = spdiags(grad_a{1}(Xq) ./ a(Xq), 0, size(Wq, 1), size(Wq, 1));
% grad_Ax{2} = spdiags(grad_a{2}(Xq) ./ a(Xq), 0, size(Wq, 1), size(Wq, 1));
% grad_Ax{3} = spdiags(grad_a{3}(Xq) ./ a(Xq), 0, size(Wq, 1), size(Wq, 1));
% 
% Mult = Nrm{1} * grad_Ax{1} + Nrm{2} * grad_Ax{2} + Nrm{3} * grad_Ax{3};
% 
% V = (1i/4) * integral(Xq, Gamma, Gxy, Vh);
% V = V - 1/2/pi * regularize(Xq, Gamma,'[log(r)]',Vh);
% 
% 
% %%%%%%%%%%%%%%%%%%% MULT %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V = Mu' * Wq * Mult * V; 
% 
% K = (1i/4) * integral(Xq, Gamma, dyGxy, ntimes(Uh));
% K = K - 1/(2*pi) * regularize(Xq, Gamma, 'grady[log(r)]',ntimes(Uh));
% 
% I = integral(Gamma, ntimes(Uh),Axm1_gradAx,Uh);
% 
% 
% %%%%%%%%%%%%%%%%%%% MULT %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K = 0.5*I - Mu' * Wq * Mult * K;
% 
% I1 = integral(Gamma, Vh, Uh);
% I2 = integral(Gamma, Vh, Vh);




end

