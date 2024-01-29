function mesh = mshCubeLevel(n, L)

% Delaunay mesh
x = [-L 0 L];
y = [-L 0 L];
z = [-L 0 L];

if n > 1

    for level = 1:(n-1)
        
        x_new = zeros(1, 2*length(x)-1);
        
        I_old = 1:2:(2*length(x) - 1);
        
        I_new = setdiff(1:(2*length(x) - 1), I_old);
        
        x_new(I_old) = x;
        
        x_new(I_new) = 0.5 * (x(1:end-1) + x(2:end));
        
        x = x_new;
        
    end


end

y = x; z = x;

[x,y,z] = meshgrid(x,y,z);
DT      = delaunayTriangulation([x(:) y(:) z(:)]);

% Build mesh
mesh = msh(DT.Points,DT.ConnectivityList);
end