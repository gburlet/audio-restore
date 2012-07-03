function [a, G] = arModel(r, p)
% Gregory Burlet
% MUMT 605
% December 9, 2011
%
% ARMODEL
% Calculates autoregressive parameters using covariance method

N = length(r);

G = toeplitz(fliplr(r'));
G = G(1:N-p,N-p+1:N);

% estimate AR coefficients and error
a = (G'*G)\G'*r(p+1:N);

end