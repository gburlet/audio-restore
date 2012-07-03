function [c, G] = sinModel(x, q)
% Gregory Burlet
% MUMT 605
% December 9, 2011
%
% This function determines the parameters for the sinusoidal basis

% get frequencies for sinusoidal basis
G = calcSinBasis(x, q);

% calculate sinusoid coefficients using least squares
c = (G'*G)\G'*x;

end

function [G] = calcSinBasis(x, q)
% return G: basis matrix of sinusoids

N = length(x);

numFreq = floor(q/2);

% N is a power of 2 for efficiency
X = fft(x, N);
  
% since real signal, we can remove second half of spectrum and DC
binEnergy = abs(X(2:N/2));
  
[~,i] = sort(binEnergy,'descend');
% get normalized frequencies
freqs = i(1:numFreq) ./ N;

% calculate basis matrix of sinusoids
numBases = 2*numFreq+1;

% TODO: create without for loop
G = zeros(N, numBases);

t = 0:N-1;
omega = 2*pi*freqs;

for k = 1:numFreq
    G(:,2*k-1) = cos(t.*omega(k));  % odd
    G(:,2*k) = sin(t.*omega(k));    % even
end

% DC
G(:,numBases) = ones(N,1);

end