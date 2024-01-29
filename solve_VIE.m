 function [uh, uext] = solve_VIE(mesh, pts, k0, kx, typ, tol, theta)



Omega = dom(mesh, 3);


Vh = fem(mesh, typ);



G = @(X, Y) femGreenKernel(X, Y, '[H0(kr)]', k0);

uinc = @(X) exp(1i*k0*(cos(theta)*X(:,1) +sin(theta)*X(:,2)));


if tol < 1e-12

K = (1i/4) * integral(Omega, Omega,  Vh, G, Vh);

else
    
K = (1i/4) * integral(Omega, Omega,  Vh, G, Vh, tol);
end

Kreg = -(1/(2 * pi)) * regularize(Omega,Omega,Vh,'[log(r)]',Vh);

I = integral(Omega, Vh, Vh);

F = integral(Omega, Vh, uinc);


% SOLVE FOR ONE KX %
%-----------------------------------------------------------------%
if ~iscell(kx)


if strcmp(typ, 'P0')
    
    beta = @(X) kx(X).^2 - k0^2;
    BetaP0 = integral(Omega, Vh, beta);
    Mass = integral(Omega, Vh, Vh);
    
    Beta = Mass \ BetaP0;
    
    dim = length(Beta);
    
    Beta = spdiags(Beta, 0, dim, dim); 
    
elseif strcmp(typ, 'P1')
    
    dim = size(mesh.vtx, 1);
    
    
    beta = @(X) kx(X).^2 - k0^2;
    Beta = spdiags(beta(mesh.vtx), 0, dim, dim); 
    
    if tol > 1e-12
    Beta = hmx(double(mesh.vtx),double(mesh.vtx),Beta,tol);
    end
    
end

if strcmp(typ, 'P0')
    A = (K + Kreg) * Beta;

elseif strcmp(typ, 'P1')
    
    A = (K + Kreg) * Beta;

end

M = I-A;


uh = M\F;




% Evaluates using Representation formula

K = (1i/4) * integral(pts, Omega, G, Vh);
Kreg = -(1/(2 * pi)) * regularize(pts, Omega, '[log(r)]',Vh);

aux = uinc(pts);

if strcmp(typ, 'P0')
    uext = aux + (K + Kreg)* (Beta*uh);

elseif strcmp(typ, 'P1')
    uext = aux + (K + Kreg)* (Beta*uh);
    
end


% SOLVE FOR MANY KX %
%-------------------------------------------------------------------%
else
    
Nk = length(kx);

for count = 1:Nk
    
if strcmp(typ, 'P0')
    
    beta = @(X) kx{count}(X).^2 - k0^2;
    BetaP0 = integral(Omega, Vh, beta);
    Mass = integral(Omega, Vh, Vh);
    
    Beta = Mass \ BetaP0;
    
    dim = length(Beta);
    
    Beta = spdiags(Beta, 0, dim, dim); 
    
elseif strcmp(typ, 'P1')
    
    dim = size(mesh.vtx, 1);
    
    
    beta = @(X) kx{count}(X).^2 - k0^2;
    Beta = spdiags(beta(mesh.vtx), 0, dim, dim); 
    
    if tol > 1e-12
    Beta = hmx(double(mesh.vtx),double(mesh.vtx),Beta,tol);
    end
    
end

if strcmp(typ, 'P0')
    A = (K + Kreg) * Beta;

elseif strcmp(typ, 'P1')
    
    A = (K + Kreg) * Beta;

end

M = I-A;

uh = M\F;


% Evaluates using Representation formula

KK = (1i/4) * integral(pts, Omega, G, Vh);
KKreg = -(1/(2 * pi)) * regularize(pts, Omega, '[log(r)]',Vh);

aux = uinc(pts);

if strcmp(typ, 'P0')
    uex = aux + (KK + KKreg)* (Beta*uh);

elseif strcmp(typ, 'P1')
    uex = aux + (KK + KKreg)* (Beta*uh);
    
end

uext(count) = uex;

end

end
end
