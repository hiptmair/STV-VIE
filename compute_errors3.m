%%%%%%%%%%%
% For 3D meshes it assumes that mesh.col is used appropriately:
% i.e. m.col = (1:Nelt).' and mref is refined from m using
% mshMidpointTetra.
%%%%%%%%%%%

function err = compute_errors3(m, mref, Uh, Uex, type)

if size(m, 1) ~= size(mref, 1)

[Uhh, ~] = interpolation(m, mref, Uh);

else
    
Uhh = Uh;
end

if size(m.elt, 2) == 3
domain = dom(mref, 3);

elseif size(m.elt, 2) == 4
domain = dom(mref, 4);
    
end
fe = fem(mref, 'P1');

% Mesh integration data
[X,Wx] = domain.qud;

% Finite element matrix
Mu = fe.uqm(domain);

    
Uapp = Mu * Uhh;
Uexact = Mu * Uex;

switch type
    case 'L2'
        err = sqrt(sum(Wx.*(real(Uexact- Uapp)).^2));
        
        err = real(err);
        err = err / (sqrt(sum(Wx.*(real(Uexact)).^2)));
    case 'H1'
        
        gradef = grad(fe);
        Gu = gradef.uqm(domain);
        
        if size(m.elt, 2) == 3
        DxUapp = Gu{1} * Uhh;
        DyUapp = Gu{2} * Uhh;
        
        DxUexact = Gu{1} * Uex;
        DyUexact = Gu{2} * Uex;
        
        err = sqrt(sum(Wx.*(real(Uexact-Uapp).^2+real(DxUexact-DxUapp).^2+real(DyUexact-DyUapp).^2)));
        
        err = real(err);
        err = err / real(sqrt(sum(Wx.*(real(Uexact).^2+real(DxUexact).^2+real(DyUexact).^2))));
        elseif size(m.elt, 2) == 4
           
        DxUapp = Gu{1} * Uhh;
        DyUapp = Gu{2} * Uhh;
        DzUapp = Gu{3} * Uhh;
        
        DxUexact = Gu{1} * Uex;
        DyUexact = Gu{2} * Uex;
        DzUexact = Gu{3} * Uex;
        
        err = sqrt(sum(Wx.*(real(Uexact-Uapp).^2+real(DxUexact-DxUapp).^2+real(DyUexact-DyUapp).^2+real(DzUexact-DzUapp).^2)));
        
        err = real(err);
        err = err / real(sqrt(sum(Wx.*(real(Uexact).^2+real(DxUexact).^2+real(DyUexact).^2+real(DzUexact).^2)))); 
        end
    otherwise
        error('Unknown error type. Known types are ''L2'' or ''H1''.');
end


end