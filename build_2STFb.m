function [A0, A1] = build_2STFb(mesh, k0, k, a, alpha, tol)


Gamma = dom(mesh, 3);
Uh = fem(mesh, 'P1');
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

Mu = femLagrangePn(Uh, Gamma);
[~,Pu] = Uh.unk; 

Mu = Mu * Pu;  

[Xq, Wq, elt2qud] = Gamma.qud;

Wq = spdiags(Wq, 0, size(Wq, 1), size(Wq, 1));

Ax = spdiags(a(Xq), 0, size(Wq, 1), size(Wq, 1));
mAx = spdiags(1./a(mesh.vtx), 0, size(mesh.vtx, 1), size(mesh.vtx, 1));


V = (1i/4) * integral(Xq, Gamma, Gxy, Uh);
V = V - 1/2/pi * regularize(Xq, Gamma,'[log(r)]',Uh);

V = Mu' * Wq * Ax * V; 


K = (1i/4) * integral(Xq, Gamma, dyGxy, ntimes(Uh));
K = K - 1/(2*pi) * regularize(Xq, Gamma, 'grady[log(r)]',ntimes(Uh));

K = Mu' * Wq * Ax * K * mAx; 


Kp = (1i/4) * integral(Gamma, Gamma, ntimes(Uh), dxGxy, Uh);
Kp = Kp - 1/(2*pi) * regularize(Gamma, Gamma, Uh,'grady[log(r)]',ntimes(Uh)).';

W = 1i/4 * (k0^2 * integral(Gamma, Gamma, ntimes(Uh),Gxy,ntimes(Uh)) ...
    - integral(Gamma, Gamma, nxgrad(Uh), Gxy, nxgrad(Uh)));
W = W -1/(2*pi) * (k0^2 * regularize(Gamma, Gamma, ntimes(Uh),'[log(r)]',ntimes(Uh)) ...
    - regularize(Gamma, Gamma, nxgrad(Uh),'[log(r)]',nxgrad(Uh)));

W = W * mAx;

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

V = (1i/4) * integral(Gamma, Gamma, Uh, Gxy, Uh);
V = V - 1/2/pi * regularize(Gamma, Gamma, Uh,'[log(r)]',Uh);

K = (1i/4) * integral(Gamma, Gamma, Uh, dyGxy, ntimes(Uh));
K = K - 1/(2*pi) * regularize(Gamma, Gamma, Uh,'grady[log(r)]',ntimes(Uh));

Kp = (1i/4) * integral(Gamma, Gamma, ntimes(Uh), dxGxy, Uh);
Kp = Kp - 1/(2*pi) * regularize(Gamma, Gamma, Uh,'grady[log(r)]',ntimes(Uh)).';

W = 1i/4 * (k^2 * integral(Gamma, Gamma, ntimes(Uh),Gxy,ntimes(Uh)) ...
    - integral(Gamma, Gamma, nxgrad(Uh), Gxy, nxgrad(Uh)));
W = W -1/(2*pi) * (k^2 * regularize(Gamma, Gamma, ntimes(Uh),'[log(r)]',ntimes(Uh)) ...
    - regularize(Gamma, Gamma, nxgrad(Uh),'[log(r)]',nxgrad(Uh)));
% toc

A1 = [-K V;-W Kp];

end

