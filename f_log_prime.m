function y = f_log_prime(x)
% Computes the derivative of the inverse logit link function

y = 1.0./(2 + exp(-x)+exp(x));
end