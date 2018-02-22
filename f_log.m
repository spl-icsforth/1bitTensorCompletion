function y = f_log(x)
% Computes the inverse logit link function

y = 1.0./(1+exp(-x));
end