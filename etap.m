function y = etap(x)

y = eta(x) .* (-2 ./ (1 - x.^2).^2);
I = isinf(y);

y(I) = 0;

I = (abs(1-x.^2) <= 1e-12);

y(I) = 0;

end