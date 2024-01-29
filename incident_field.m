function [uinc, gradxUinc] = incident_field(k0, k, ax, theta, type)


if strcmp(type, 'PW')
   
    uinc = @(X) exp(1i*k0* (cos(theta) * X(:, 1) + sin(theta) * X(:, 2)));
    
    gradxUinc{1} = @(X) 1i*k0*cos(theta) .* uinc(X);
    gradxUinc{2} = @(X) 1i*k0*sin(theta) .* uinc(X);
    gradxUinc{3} = @(X) 1i*k0*0 .* uinc(X);

    
elseif strcmp(type, 'SRC')

    X0 = [3/2 3/2];
    uinc = @(X) 1i / 4 * besselh(0, k0 * sqrt( (X(:, 1) - X0(1)).^2 + (X(:, 2) - X0(2)).^2 ));
    
    gradxUinc{1} = @(X) -1i / 4 * k0 * besselh(1, k0 * sqrt( (X(:, 1) + X0(1)).^2 + (X(:, 2) + X0(2)).^2)) .* (X(:, 1) + X0(1)) ./ sqrt((X(:, 1) + X0(1)).^2 + (X(:, 2) + X0(2)).^2);
    gradxUinc{2} = @(X) -1i / 4 * k0 * besselh(1, k0 * sqrt( (X(:, 1) + X0(1)).^2 + (X(:, 2) + X0(2)).^2)) .* (X(:, 2) + X0(2)) ./ sqrt((X(:, 1) + X0(1)).^2 + (X(:, 2) + X0(2)).^2);
    gradxUinc{3} = @(X) 0 .* X(:, 1);

    
    
elseif strcmp(type, 'Exact')
    

    Theta = @(X) mod(atan2(X(:, 2), X(:, 1)), 2*pi);
    Rad = @(X) X(:, 1).^2 + X(:, 2).^2;

    dRad{1} = @(X) 2 * X(:, 1);
    dRad{2} = @(X) 2 * X(:, 2);
    dRad{3} = @(X) 0 * X(:, 3);

    dTheta{1} = @(X) -sin(Theta(X)) ./ Rad(X).^(1/2);
    dTheta{2} = @(X)  cos(Theta(X)) ./ Rad(X).^(1/2);
    dTheta{3} = @(X) 0*X(:, 1);


    J = @(X) k.^(-1/2) .* besselj(2 / 3, k * sqrt(Rad(X)));

    dJ{1} = @(X) k^(1/2) * 0.5 * (besselj(-1/3, k * sqrt(Rad(X))) - ...
                                  besselj(5/3, k * sqrt(Rad(X))) ) .* ...
                                  dRad{1}(X) ./ sqrt(Rad(X))/ 2;



    dJ{2} = @(X) k^(1/2) * 0.5 * (besselj(-1/3, k * sqrt(Rad(X))) - ...
                                  besselj(5/3, k * sqrt(Rad(X))) ) .* ...
                                  dRad{2}(X) ./ sqrt(Rad(X)) /2;


    Uin = @(X) J(X) .* sin(2 / 3 * Theta(X));
    Uex = @(X) 1i / 4 * besselh(0, k0 * sqrt( (X(:, 1) + 1/4).^2 + (X(:, 2) + 1/4).^2) );


    gradxUin{1} = @(X) sin(2 / 3 * Theta(X)) .* dJ{1}(X) + ...
                       J(X) .* cos(2 / 3 * Theta(X)) .* dTheta{1}(X) * 2 / 3;

    gradxUin{2} = @(X) sin(2 / 3 * Theta(X)) .* dJ{2}(X) + ...
                       J(X) .* cos(2 / 3 * Theta(X)) .* dTheta{2}(X) * 2 / 3;
    gradxUin{3} = @(X) 0 .* X(:, 1);


    gradxUex{1} = @(X) -1i / 4 * k0 * besselh(1, k0 * sqrt( (X(:, 1) + 1/4).^2 + (X(:, 2) + 1/4).^2)) .* (X(:, 1) + 1/4) ./ sqrt((X(:, 1) + 1/4).^2 + (X(:, 2) + 1/4).^2);
    gradxUex{2} = @(X) -1i / 4 * k0 * besselh(1, k0 * sqrt( (X(:, 1) + 1/4).^2 + (X(:, 2) + 1/4).^2)) .* (X(:, 2) + 1/4) ./ sqrt((X(:, 1) + 1/4).^2 + (X(:, 2) + 1/4).^2);
    gradxUex{3} = @(X) 0 .* X(:, 1);


    uinc = @(X) Uex(X) - Uin(X);

    gradxUinc{1} = @(X) gradxUex{1}(X) - gradxUin{1}(X);
    gradxUinc{2} = @(X) gradxUex{2}(X) - gradxUin{2}(X);
    gradxUinc{3} = @(X) 0 .* X(:, 1);


end
end