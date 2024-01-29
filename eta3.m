function y = eta3(rad, X)


r2x = rad^2 * X(:, 1).^2;
r2y = rad^2 * X(:, 2).^2;

y = exp(2) * exp(-1./(1-r2x)) .* exp(-1./(1-r2y));

I = (r2x >= 1);

y(I) = 0;


I = (r2y >= 1);

y(I) = 0;

end
