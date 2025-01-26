function fit = getListFit(fit,list,M)
%GETLISTFIT Summary of this function goes here
%   Detailed explanation goes here
fit.PeakList = list.peakList;
fit.CompNams = string(list.compNam);

iList = round(list.peakList)==M;
n = sum(iList);

ind = fit.param.massRange==M;
iM = round(fit.mz)==M;

fit.peaks(:,ind) = nan;
fit.peaks(1:n,ind) = fit.PeakList(iList);

list.compNams = fit.CompNams;
isoSpecs = getIsoYs(fit,list,M,true);
specs = fit.specs(iM,:)-isoSpecs;
param = rmfield(fit.param,'BL');
param.BL = fit.param.BL(ind);
param.baseline = 0;
param.peakList = fit.PeakList(iList);
[~, ~, ~, H, ~, bl] = fitPeaks(fit.mz(iM), specs, [], param);

fit.H(1:n,:,ind) = H;
if n < param.peakRange(end)
    fit.H(n+1:end,:,ind) = nan;
end
fit.specsFinal(iM,:) = specs;
fit.param.BL(ind) = bl;
fit.isoSpecs(iM,:) = isoSpecs;
end

