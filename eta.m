function y = eta(x)

y = exp(1) * exp(-1./(1-x.^2));
I = isinf(y);

y(I) = 0;

I = (x.^2 > 1);

y(I) = 0;

end