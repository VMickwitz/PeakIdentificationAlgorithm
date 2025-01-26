function [prefit] = getPrelimFit(fit,M)
%GETPRELIMFIT Returns a fit structure for the preliminary fit at mass M.

fitpre = fit.fitPrelim;

ind = fit.param.massRange == M;

Apre = fit.param.A(ind);
%fitpre.Alims(:,2,ind)
ia = find(fitpre.Alims(:,2,ind) > Apre,1,'first');

fitpre.peaks = fitpre.peaks(ia,:,ind);
fitpre.H = fitpre.H(:,:,ia,ind);
fitpre.numbers = fitpre.numbers(ia,ind);
fitpre.param = fit.param;
fitpre.param.massRange = M;
fitpre.param.W = fit.param.W(ind);
fitpre.param.BL = fitpre.BLs(ia,ind)*fitpre.param.W;
fields = ["Avec","removed","chi","BLs","Alims","baseline",...
    "nPeaks","chis"];
prefit = rmfield(fitpre,fields);

prefit.mz = fit.mz;
prefit.specs = fit.specs;
prefit.isoSpecs = nan(size(prefit.specs));
end

