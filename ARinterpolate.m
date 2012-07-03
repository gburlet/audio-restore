function [xi] = ARinterpolate(x, a, i)
% Gregory Burlet
% MUMT 605
% December 9, 2011
%
% This function estimates missing data using the autoregressive model

numMissing = sum(i == 1);
xi = zeros(numMissing,1);

if numMissing > 0
    N = length(x);
    p = length(a);
    
    A = triu(toeplitz([fliplr(-a'), 1, zeros(1,N-p-1)]));
    A = A(1:N-p,:);
    
    % get known values
    x_i = x(i == 0);
    
    Ai = A(:, i == 1);
    A_i = A(:, i == 0);
    
    xi = -(Ai'*Ai)\Ai'*A_i*x_i;
end

end

