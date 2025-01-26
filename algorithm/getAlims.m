function [lims,peaks,h_outs,numbers,chis,removed,BLs] = getAlims(x,y,param)
% Returns the limits of A within which a fit is obtained, as well as the
% corresponding fits.
%   outputs:
%       lims:       peakn*2 array with the lower and upper limits of A.
%       peaks:      The peaks that the corresponding values of A give.
%       h_outs:     Heights that the corresponding values of A give.
%       numbers:    Number of peaks fit at corresponding A.
%       chis:       Chi values of all peak numbers
%       removed:    Spectra removed.
%       BLs:        Fitted baselines for corresponding values of A.

param.method = 3;
param.returnAll = true;
if isfield(param,'initA')
    param.A = param.initA;
else
    param.A = 1;
end
param.useA = false;
param.findA = false;
if isfield(param,'Rtol')
    param = rmfield(param,'Rtol');
end

[numbers, score, peaksfit, removed, h_outs, BLs] = guessPeakNumber(x, y, param);

[numbers,i_sort] = sort(numbers);
score = score(i_sort);

% Calculate chi values from the score.
chis = (score - numbers*log(length(x)*size(y,2)))/param.A;
chi_out = chis(:);


[Alims,numbers,~,irm] = getLims(numbers,chis);

Alims = Alims.*log(length(x)*size(y,2));

BLs(irm,:) = [];
peaksfit(irm,:) = [];
h_outs(:,:,irm) = [];
chi_out(irm,:) = [];

chis = chi_out;
lims = Alims;
peaks = peaksfit;

end

