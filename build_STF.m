function A = build_STF(mesh, k0, k, tol)


Gamma = dom(mesh, 3);
Vh = fem(mesh, 'P0');
Uh = fem(mesh, 'P1');
%%


Gxy = @(X,Y) femGreenKernel(X, Y, '[H0(kr)]', k0);
dyGxy{1} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]1', k0);
dyGxy{2} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]2', k0);
dyGxy{3} = @(X,Y) femGreenKernel(X, Y, 'grady[H0(kr)]3', k0);

dxGxy{1} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]1', k0);
dxGxy{2} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]2', k0);
dxGxy{3} = @(X,Y) femGreenKernel(X, Y, 'gradx[H0(kr)]3', k0);


if tol > 1e-12
%%
% disp('--------------------------------')
% disp('Assembly of Integral Operator A0');
% tic

V = (1i/4) * integral(Gamma, Gamma, Vh, Gxy, Vh, tol);
V = V - 1/2/pi * regularize(Gamma, Gamma, Vh,'[log(r)]',Vh);

K = (1i/4) * integral(Gamma, Gamma, Vh, dyGxy, ntimes(Uh), tol);
K = K - 1/(2*pi) * regularize(Gamma, Gamma, Vh,'grady[log(r)]',ntimes(Uh));

Kp = (1i/4) * integral(Gamma, Gamma, ntimes(Uh), dxGxy, Vh, tol);
Kp = Kp - 1/(2*pi) * regularize(Gamma, Gamma, Vh,'grady[log(r)]',ntimes(Uh)).';

W = 1i/4 * (k0^2 * integral(Gamma, Gamma, ntimes(Uh),Gxy,ntimes(Uh), tol) ...
    - integral(Gamma, Gamma, nxgrad(Uh), Gxy, nxgrad(Uh), tol));
W = W -1/(2*pi) * (k0^2 * regularize(Gamma, Gamma, ntimes(Uh),'[log(r)]',ntimes(Uh)) ...
    - regularize(Gamma, Gamma, nxgrad(Uh),'[log(r)]',nxgrad(Uh)));
% toc


% A0 = bmm({-W, Kp; -K, V});
% A0 = bmm({-K, V; -W, Kp});
% A0 = [-K V;-W Kp];
A0 = [-W Kp; -K V];

 
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

V = (1i/4) * integral(Gamma, Gamma, Vh, Gxy, Vh, tol);
V = V - 1/2/pi * regularize(Gamma, Gamma, Vh,'[log(r)]',Vh);

K = (1i/4) * integral(Gamma, Gamma, Vh, dyGxy, ntimes(Uh), tol);
K = K - 1/(2*pi) * regularize(Gamma, Gamma, Vh,'grady[log(r)]',ntimes(Uh));

Kp = (1i/4) * integral(Gamma, Gamma, ntimes(Uh), dxGxy, Vh, tol);
Kp = Kp - 1/(2*pi) * regularize(Gamma, Gamma, Vh,'grady[log(r)]',ntimes(Uh)).';

W = 1i/4 * (k^2 * integral(Gamma, Gamma, ntimes(Uh),Gxy,ntimes(Uh), tol) ...
    - integral(Gamma, Gamma, nxgrad(Uh), Gxy, nxgrad(Uh), tol));
W = W -1/(2*pi) * (k^2 * regularize(Gamma, Gamma, ntimes(Uh),'[log(r)]',ntimes(Uh)) ...
    - regularize(Gamma, Gamma, nxgrad(Uh),'[log(r)]',nxgrad(Uh)));
% toc


% A1 = bmm({-W, Kp; -K, V});
% A1 = bmm({-K, V; -W, Kp});
% A1 = [-K V;-W Kp];
A1 = [-W Kp; -K V];

A = A0+A1;

else

% disp('--------------------------------')
% disp('Assembly of Integral Operator A0');
% tic

V = (1i/4) * integral(Gamma, Gamma, Vh, Gxy, Vh);
V = V - 1/2/pi * regularize(Gamma, Gamma, Vh,'[log(r)]',Vh);

K = (1i/4) * integral(Gamma, Gamma, Vh, dyGxy, ntimes(Uh));
K = K - 1/(2*pi) * regularize(Gamma, Gamma, Vh,'grady[log(r)]',ntimes(Uh));

Kp = (1i/4) * integral(Gamma, Gamma, ntimes(Uh), dxGxy, Vh);
Kp = Kp - 1/(2*pi) * regularize(Gamma, Gamma, Vh,'grady[log(r)]',ntimes(Uh)).';

W = 1i/4 * (k0^2 * integral(Gamma, Gamma, ntimes(Uh),Gxy,ntimes(Uh)) ...
    - integral(Gamma, Gamma, nxgrad(Uh), Gxy, nxgrad(Uh)));
W = W -1/(2*pi) * (k0^2 * regularize(Gamma, Gamma, ntimes(Uh),'[log(r)]',ntimes(Uh)) ...
    - regularize(Gamma, Gamma, nxgrad(Uh),'[log(r)]',nxgrad(Uh)));
% toc


% A0 = bmm({-W, Kp; -K, V});
% A0 = bmm({-K, V; -W, Kp});
% A0 = [-K V;-W Kp];
A0 = [-W Kp; -K V];

 
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


% A1 = bmm({-W, Kp; -K, V});
% A1 = bmm({-K, V; -W, Kp});
% A1 = [-K V;-W Kp];
A1 = [-W Kp; -K V];

A = A0+A1;

end

