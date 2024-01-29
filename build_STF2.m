function [A0, A1] = build_STF2(mesh, k0, k, a, alpha, tol)


Gamma = dom(mesh, 3);
Vh = fem(mesh, 'P0');
Uh = fem(mesh, 'P1');
nUh = ntimes(Uh);
curlUh = nxgrad(Uh);
%%


Gxy = @(X,Y) femGreenKernel(X, Y, '[H0(kr)]', k0);

comm_Gxy = @(X,Y) (a(Y) - a(X)) .* femGreenKernel(X, Y, '[H0(kr)]', k0);

dyGxy{1} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]1', k0);
dyGxy{2} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]2', k0);
dyGxy{3} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]3', k0);


comm_dyGxy{1} = @(X,Y) (a(Y) - a(X)) .* femGreenKernel(X, Y, 'grady[H0(kr)]1', k0);
comm_dyGxy{2} = @(X,Y) (a(Y) - a(X)) .* femGreenKernel(X, Y, 'grady[H0(kr)]2', k0);
comm_dyGxy{3} = @(X,Y) (a(Y) - a(X)) .* femGreenKernel(X, Y, 'grady[H0(kr)]3', k0);


dxGxy{1} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]1', k0);
dxGxy{2} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]2', k0);
dxGxy{3} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]3', k0);


comm_dxGxy{1} = @(X,Y) 1./ a(X) .* (a(Y) - a(X)) .* femGreenKernel(X, Y, 'gradx[H0(kr)]1', k0);
comm_dxGxy{2} = @(X,Y) 1./ a(X) .* (a(Y) - a(X)) .* femGreenKernel(X, Y, 'gradx[H0(kr)]2', k0);
comm_dxGxy{3} = @(X,Y) 1./ a(X) .* (a(Y) - a(X)) .* femGreenKernel(X, Y, 'gradx[H0(kr)]3', k0);



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
mAx = spdiags(1./a(Xq), 0, size(Wq, 1), size(Wq, 1));


V = (1i/4) * integral(Xq, Gamma, Gxy, Vh);
V = V - 1/2/pi * regularize(Xq, Gamma,'[log(r)]',Vh);

V = Mv' * Wq * Ax * V; 

V = V + (1i/4) * integral(Gamma, Gamma, Vh, comm_Gxy, Vh);


K = (1i/4) * integral(Gamma, Gamma, Vh, dyGxy, ntimes(Uh));
K = K - 1/(2*pi) * regularize(Gamma, Gamma, Vh,'grady[log(r)]',ntimes(Uh));

Kp = (1i/4) * integral(Gamma, Gamma, ntimes(Uh), dxGxy, Vh);
Kp = Kp - 1/(2*pi) * regularize(Gamma, Gamma, Vh,'grady[log(r)]',ntimes(Uh)).';

Kp = Kp + (1i/4) * integral(Gamma, Gamma, ntimes(Uh), comm_dxGxy, Vh);


W1 = integral(Xq, Gamma, Gxy, ntimes(Uh));
W2 = integral(Xq, Gamma, Gxy, nxgrad(Uh));

reg1 = regularize(Xq, Gamma,'[log(r)]',ntimes(Uh));
reg2 = regularize(Xq, Gamma,'[log(r)]',nxgrad(Uh));

W1{1} = 1i/4 * k0^2 * W1{1} -1/(2*pi) * k0^2 * reg1{1};
W1{2} = 1i/4 * k0^2 * W1{2} -1/(2*pi) * k0^2 * reg1{2};
W1{3} = 1i/4 * k0^2 * W1{3} -1/(2*pi) * k0^2 * reg1{3};

W2{1} =-1i/4 * W2{1} +1/(2*pi) * reg2{1};
W2{2} =-1i/4 * W2{2} +1/(2*pi) * reg2{2};
W2{3} =-1i/4 * W2{3} +1/(2*pi) * reg2{3};


W1 = Mnu{1}' * Wq * mAx * W1{1} + ...
     Mnu{2}' * Wq * mAx * W1{2}+ ...
     Mnu{3}' * Wq * mAx * W1{3};


W2 = Mcu{1}' * Wq * mAx * W2{1} + ...
     Mcu{2}' * Wq * mAx * W2{2}+ ...
     Mcu{3}' * Wq * mAx * W2{3};

W = W1 + W2;

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
alphaI = integral(Gamma, Uh, alpha, Vh);
A1 = [-K V;-W Kp-alphaI];

end

