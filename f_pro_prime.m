function y = f_pro_prime(x, mu, sigma)
% Computes the pdf of a Gaussian distribution at specified values

y = exp(-0.5 .* ((x - mu) ./ sigma).^2) ./ (sqrt(2 * pi) .* sigma);
end