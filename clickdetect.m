function [i] = clickdetect(ehat, stdhat, detThresh, detStretch)
% Gregory Burlet
% MUMT 605
% December 9, 2011
%
% returns indices of detected noise

% widen so indices don't go out of bounds.
i_temp = [zeros(detStretch,1); abs(ehat) > detThresh*stdhat; zeros(detStretch,1)];
i = zeros(length(ehat),1);

% bad, bad for loop
% in future find way to vectorize this widening process
for k = 1:length(ehat)
    i(k) = any(i_temp(k:k+2*detStretch));
end

end