function y = eta2(rad, X)

r2 = rad^2 * (X(:, 1).^2 + X(:, 2).^2 + X(:, 3).^2);

y = exp(1) * exp(-1./(1-r2));

I = (r2 >= 1);

y(I) = 0;

end

