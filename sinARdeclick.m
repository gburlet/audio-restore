function [x, i] = sinARdeclick(x, p, q, N, H, detThresh, detStretch, numIter)
% Gregory Burlet
% MUMT 605
% December 9, 2011
%
% This function removes clicks from the specified wavfile using a
% sinusoid basis and stochastic autoregressive residual model.
%
% param x: matrix of corrupted audio data (columns = channels, rows = samples)
% param p: AR model order
%          default = 31
% param q: number of sinusoids in basis
%          default = 31
% param N: length of window/frame (samples)
%          default = 2048
% param H: window overlap (samples)
%          default = N
% param detThresh: threshold for click detection (number
%          default = 4
% param detStretch: amount to stretch detection boundaries
%          default = 4
% param numIter: amount of iterations for iterative optimizations
%          default = 3
%
% return x: declicked audio data
% return i: binary click indicator vector corresponding to x

% set default values
if nargin < 8
    numIter = 3;
end
if nargin < 7
    detStretch = 4;
end
if nargin < 6
    detThresh = 4;
end
if nargin < 4
    N = 2048;
end
if nargin < 5
    H = floor(N/2);
end
if nargin < 3
    q = 31;
end
if nargin < 2
    p = 31;
end
if nargin == 0
    error('Please load in a sound file');
end

% make window a power of two for efficiency
N = 2^nextpow2(N)-p;

x_len = size(x,1);
numWin = floor(x_len/H);

numChannels = size(x,2);

% perform declicking on each channel seperately
for chan = 1:numChannels
    % init indicator vector
    i = zeros(x_len,numChannels);
    
    % for each frame of signal
    for k = 0:numWin
        % shift over by p to accompany xhead
        f_start = k*H+1+p;
        f_end = k*H+N+p;
        
        % zero pad frame if necessary
        if f_end > x_len
            frame = [x(f_start:x_len, chan); zeros(N-x_len+f_start-1,1)];
        else
            frame = x(f_start:f_end, chan);
        end
        
        % prepend last p samples from previous block
        fhead = x(f_start-p:f_start-1, chan);
        
        [clean_frame, frame_i] =  sinARdeclickFrame(frame, fhead, p, q, detThresh, detStretch, numIter);
        
        display(['status: %', num2str(f_end/(H*numWin+N+p)*100)]);
        
        % add frame_i to global indicator vector
        i(f_start:f_end, chan) = frame_i;
        
        % replace corrupted data with clean interpolated data
        x(f_start:f_end, chan) = clean_frame;
        
    end
end

end

function [x, i] =  sinARdeclickFrame(x, xhead, p, q, detThresh, detStretch, numIter)
% Gregory Burlet
% MUMT 605
% December 9, 2011
%
% Performs iterative click detection and removal on a single frame using a
% sinusoid basis and a stochastic autoregressive residual model.

N = length(x);

xall = [xhead; x];

[c, Gsin] = sinModel(xall, q);

% calculate residuals
sins = Gsin*c;
% r = x - Gc
r = xall - sins;

% calculate AR parameters from residuals
[a, Gar] = arModel(r, p);

% innovations e = (x_1 - Ga)
ehat = r(p+1:N+p) - Gar*a;

% calculate robust estimate of std deviation
stdhat = median(abs(ehat))./0.6745;

% estimate click locations
i = clickdetect(ehat, stdhat, p, detThresh, detStretch);
% accompany head of length p
iall = [zeros(p,1); i];
    
if sum(i == 1) > 0
    % iterative optimization
    for k = 1:numIter
        % interpolate unknown values of x: x(i)
        ri = ARinterpolate(r, a, iall);
        
        % add back in sinusoids to interpolated residual data
        % update frame with interpolated values
        x(i == 1) = ri + sins(iall == 1);
        xall = [xhead; x];
        
        % now estimate sinusoid coefficients again
        [c, Gsin] = sinModel(xall, q);
        
        sins = Gsin*c;
        r = xall - sins;
        
        % now estimate AR parameters again
        [a, ~] = arModel(r, p);
    end
    
end

end