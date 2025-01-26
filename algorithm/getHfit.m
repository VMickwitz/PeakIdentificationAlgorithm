function [H_out,bl,dmc,Rfit] = getHfit(fit,peaks,M)
% Fits Heights and baselines to fit structure at mass M

if isfield(fit.param,'BL')
    par = rmfield(fit.param,'BL');
else
    par = fit.param;
end
if isfield(fit.param,'dmfit')
    par = rmfield(par,'dmfit');
end

par.peakList = peaks(round(peaks)==M);

ind = round(fit.mz) == M;
[~,~,~,H_out,Rfit,bl,dmc] = fitPeaks(fit.mz(ind), fit.specs(ind,:), [], par);

end

