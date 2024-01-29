function [vec2, coefs] = interpolation(mesh, mesh2, vec1)

if size(mesh.elt, 2) == 3

vec2 = zeros(size(mesh2.vtx, 1), 1);

fe = fem(mesh, 'P1');

[Xdof,elt2dof] = fe.dof;

Nelts = size(elt2dof, 1);

for j = 1:Nelts
    
    idx = elt2dof(j, :);
    
    elt = Xdof(idx, :);
    
    vals = vec1(idx);
    
    elt(4, :) = elt(1, :);
    
  
   
    in = inpolygon(mesh2.vtx(:, 1), mesh2.vtx(:, 2), elt(:, 1), elt(:, 2));
    
    Y = mesh2.vtx(in, 1:2);
    
%     plot(elt(:, 1), elt(:, 2), 'LineWidth', 1.5);hold on;
%     scatter(Y(:, 1), Y(:,2), 'filled'); 
%     title(num2str(j));hold off;
%     
    elt = elt(:, 1:2);
    
    rhs = (Y - elt(1, 1:2))';
    
    E1 = elt(2, :) - elt(1, :);
    E2 = elt(3, :) - elt(1, :);
    
    J = [E1; E2]';
    
    coefs = J \ rhs;
    
    
%     pts2 = elt(1, :)' + coefs(1, :) .* E1' + coefs(2, :) .* E2';
%     
%     diff = pts2' - Y;
    
    coefs = coefs';
    
    vec2(in) = (1- coefs(:, 1) - coefs(:, 2)) * vals(1) + ...
                vals(2) * coefs(:, 1) + vals(3) * coefs(:, 2);
        
    
end

elseif size(mesh.elt, 2) == 2

vec2 = zeros(size(mesh2.vtx, 1), 1);

fe = fem(mesh, 'P1');

[Xdof,elt2dof] = fe.dof;

Nelts = size(elt2dof, 1);

for j = 1:Nelts
    
    idx = elt2dof(j, :);
    
    elt = Xdof(idx, :);
    
    vals = vec1(idx);
    

    A = mesh2.vtx - elt(1, :);
    B = elt(2, :) - mesh2.vtx;
    E = elt(2, :) - elt(1, :);
    
    normA = sqrt(A(:, 1).^2 + A(:, 2).^2);
    
    normB = sqrt(B(:, 1).^2 + B(:, 2).^2);

    normE = norm(E);
    
    in = (normA <= normE) .* (normB <= normE);
    in = ~~in;
    
    A = A(in, :);
    
    coefs = dot(A, ones(sum(in), 1) * E, 2) / normE^2;
    
    
    
      
    vec2(in) = (1- coefs(:, 1)) * vals(1) + vals(2) * coefs(:, 1);
        
    
end

elseif size(mesh.elt, 2) == 4


vec2 = zeros(size(mesh2.vtx, 1), 1);

fe = fem(mesh, 'P1');

[Xdof,elt2dof] = fe.dof;

Nelts = size(elt2dof, 1);


X1 = mesh.vtx;

X2 = mesh2.vtx;

[~,I1,I2] = intersect(X1, X2,'rows','stable');
% P = sparse(I1,I2,1,length(X1),length(X2));

[edg, elt2edg] = mesh.edg;

X1edg = edg.ctr;

[~,I1edg,I2edg] = intersect(X1edg, X2,'rows','stable');
% P2 = sparse(I1edg,I2edg,1,length(X1edg),length(X2));

% col2 = mesh2.col;

for j = 1:Nelts
    
    idx = elt2dof(j, :);
    
    elt = Xdof(idx, :);
    
    vals = vec1(idx);
    
    
    edges = elt2edg(j, :);
    
    in = unique([I2(idx); I2edg(edges)]);
    
    Y = mesh2.vtx(in, :);
    
    rhs = (Y - elt(1, :))';
    
    E1 = elt(2, :) - elt(1, :);
    E2 = elt(3, :) - elt(1, :);
    E3 = elt(4, :) - elt(1, :);
    
    J = [E1; E2; E3]';
    
    coefs = J \ rhs;
    
    
%     pts2 = elt(1, :)' + coefs(1, :) .* E1' + coefs(2, :) .* E2' + coefs(3, :) .* E3';
%     
%     diff = pts2' - Y;
    
    coefs = coefs';
    
    vec2(in) = (1- coefs(:, 1) - coefs(:, 2) - coefs(:, 3)) * vals(1) + ...
                vals(2) * coefs(:, 1) + vals(3) * coefs(:, 2) + vals(4) * coefs(:, 3);
        
    
end


end
