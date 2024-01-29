function y = etap2(rad, X, jj)

r2 = rad^2 * (X(:, 1).^2 + X(:, 2).^2 + X(:, 3).^2);

y = eta2(rad, X) .* (-2 * rad * X(:, jj)) ./ ((1 - r2).^2);


I = (r2 >= 1);

y(I) = 0;

end