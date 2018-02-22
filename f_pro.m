function p = f_pro(x,mu,sigma)
% Computes the cdf of a Gaussian distribution at specified values 
% Inputs:  x - array of values
%          mu - mean of Gaussian
%          sigma - standard deviation of Gaussian
% Outputs: p - values of the CDF at x

% Standardize variables
z = (x-mu)./sigma;

% Compute values of CDF
p = 0.5 .* erfc(-z./sqrt(2));
end