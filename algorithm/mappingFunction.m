function [A_prime] = mappingFunction(A,W,inverse)
% maps A to more comparable A_prime. Inverts if inverse.
% Aprime is more comparable between peaks with signal area of different
% orders of magnitude. However, the case may be different with uncertain
% peakshape or mass calibration, etc. This is a key area for improvement.

% Currently it seems like high signal peaks often get overfit.
% This may be due to the following reasons:
%   1. Peak shape inaccurate.
%   2. Peak shape deformed by data aquisition process.
%   3. Mass calibration off.
%   4. Algorithm overestimates the gain in information from higher signal.
%
% Essentially anything which results in an uncertainty in signal height not
% proportional to sqrt(signal) will have an impact.
%
% Adding some relative uncertainty to the peakShape, and including it into
% calculation of chi values in fitPeaks, may improve this significantly.
% I consider it a high priority improvement.

b = 2.7;
%Wfun = @(W) W.*exp(-b*(1 + exp(-W.^0.5)));
Wfun = @(W) W;
%Wfun = @(W) 1; % Mapping disabled.

if nargin < 3 || ~inverse
    A_prime = A.*Wfun(W);
else
    A_prime = A./Wfun(W);
end

end

