function y = etap3(rad, X, jj)


r2 = rad^2 * X(:, jj).^2;

y = eta3(rad, X) .* (-2 * rad * X(:, jj)) ./ ((1 - r2).^2);

I = (r2 >= 1);

y(I) = 0;



end