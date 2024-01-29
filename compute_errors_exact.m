function err = compute_errors_exact(mesh,Uh,Uex, gradUex, type)

domain = dom(mesh, 3);

fe = fem(mesh, 'P1');
% Mesh integration data
[X,Wx] = domain.qud;

% Finite element matrix
Mu = fe.uqm(domain);

% Function to be applied
Uexact = Uex(X);
% Uexact = Uex;
if (size(Uexact,2) > 1)
    error('domDifference.m : unavailable case')
end
    
Uapp = Mu * Uh;
switch type
    case 'L2'
        err = sqrt(sum(Wx.*(real(Uexact-Uapp)).^2));
        
        err = real(err);
        err = err / (sqrt(sum(Wx.*(real(Uexact)).^2)));
    case 'H1'
        eps = 1e-6;
        gradef = grad(fe);
        Gu = gradef.uqm(domain);
        DxUapp = Gu{1} * Uh;
        DyUapp = Gu{2} * Uh;
        
%         
%         DxUexact = Gu{1} * Uex_vtx;
%         DyUexact = Gu{2} * Uex_vtx;
%         DzUexact = Gu{3} * Uex_vtx;
%         n1 = size(X,1);
%         eps1 = eps*ones(n1,1)*[1,0,0];
%         eps2 = eps*ones(n1,1)*[0,1,0];
%         eps3 = eps*ones(n1,1)*[0,0,1];
%         DxUexact = (Uex(X+eps1)-Uex(X-eps1))/(2*eps);
%         DyUexact = (Uex(X+eps2)-Uex(X-eps2))/(2*eps);
%         DzUexact = (Uex(X+eps3)-Uex(X-eps3))/(2*eps);
        
        DxUexact = gradUex{1}(X);
        DyUexact = gradUex{2}(X);
        err = sqrt(sum(Wx.*(real(Uexact-Uapp).^2+real(DxUexact-DxUapp).^2+real(DyUexact-DyUapp).^2)));
        
        err = real(err);
        err = err / real(sqrt(sum(Wx.*(real(Uexact).^2+real(DxUexact).^2+real(DyUexact).^2))));
    otherwise
        error('Unknown error type. Known types are ''L2'' or ''H1''.');
end


end