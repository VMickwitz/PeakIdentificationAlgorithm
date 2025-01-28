function [fit] = fitForA(gen,param)

% Fits peaks for the purpose of determining a suitable value for A

mz = gen.mz;
specs = gen.specs;

% Check and initialize parameters:
blvector = false;
if isfield(gen,'param')
    genParam = gen.param;
    param.H_scale = genParam.H_scale;
    param.massRange = genParam.massRange;
    massRange = param.massRange(1):param.massRange(end);
else
    if isfield(param,'massRange')
        massRange = param.massRange;
    else
        massRange = round(mz(1)):round(mz(end));
        param.massRange = massRange;
    end
end
n = length(massRange);
param.findA = false;
W = nan(n,1);
calcW = true;

if ~isfield(param,'R')
    if isfield(param,'fwhm')
        param.R = @(m) m/param.fwhm;
    else
        warning('No given resolution, using generated R')
        param.R = genParam.genR;
    end
end
if ~isfield(param,'H_scale')
    param.H_scale = 1;
end
if ~isfield(param,'baseline')
    param.baseline = genParam.noise;
elseif size(param.baseline) == size(specs)
    baseline = param.baseline;
    blvector = true;
end
if ~isfield(param,'Area')
    param.Area = true;
    % If true, peak size parameters correspond to area of peak, not height.
end
if isfield(param, 'peakShape') && ~isfield(param.peakShape, 'sp')
    pp = spline(param.peakShape.dat(:,1),param.peakShape.dat(:,2));
    param.peakShape.sp = pp;
elseif ~isfield(param, 'peakShape')
    warning("No given peakShape, using generated peakShape.")
    param.peakShape = genParam.genShape;
end
if ~isfield(param,'A_method')
    param.A_method = 1;
end
if ~isfield(param,'peakRange')
    param.peakRange = genParam.peakRange;
end

if isfield(param,'aveTime')
    specs = specs.*param.aveTime;
end

% Initialize some variables
pmin = min(param.peakRange);
pmax = max(param.peakRange);

lims = nan(pmax+1-pmin,2,n);
peaks_fit = nan(pmax+1-pmin,pmax,n);
h_fit = nan(pmax,size(specs,2),pmax,n);
numbers = nan(pmax+1-pmin,n);
removed = nan(size(specs,2),n);
chi = nan((pmax+1-pmin),n);
BLs = nan((pmax+1-pmin),n);

% Loop through all masses
h = waitbar(0,'Getting fit results...');
for i = 1:n
    
    ind = round(mz)==massRange(i);
    x = mz(ind);
    y = specs(ind,:);
    if blvector
        param.baseline = baseline(ind,:);
    end
    
    if calcW
        W(i) = mean(trapz(x,y));
        param.W = W(i);
        param.initA = param.W*exp(-2.5*(1+exp(-sqrt(param.W))));
        y = y./param.W;
        param.baseline = param.baseline./param.W;
    end
    % Make fits:
    [lims_i,peaks_fit_i,h_fit_i,numbers_i,chi_i,removed_i,BL_i] = getAlims(x,y,param);
    
    % Save results
    pn = size(numbers_i,2);
    lims(1:pn,:,i) = lims_i;
    peaks_fit(1:pn,1:pmax,i) = peaks_fit_i;
    h_fit(1:pmax,:,1:pn,i) = h_fit_i.*param.W;
    numbers(1:pn,i) = numbers_i;
    removed(1:size(removed_i,1),i) = removed_i;
    BLs(1:pn,i) = BL_i;
    %chi(1:(pmax+1-pmin),i) = chi_i';
    chi(1:pn,i) = chi_i';
    baseline(ind,:) = param.baseline;
    waitbar(i/n)
end
try
    close(h)
catch
    F = findall(0,'type','figure','tag','TMWWaitbar');
    delete(F);
end

% Add results to fit structure
fit.mz = mz;
fit.specs = specs;
fit.fitPrelim.peaks = peaks_fit;
fit.fitPrelim.H = h_fit;
fit.fitPrelim.Alims = lims;
fit.fitPrelim.numbers = numbers;
fit.fitPrelim.removed = removed;
fit.fitPrelim.chi = chi;
fit.fitPrelim.BLs = BLs;
fit.fitPrelim.baseline = baseline;
fit.param = param;
fit.param.baseline = baseline;
fit.param.W = W;
end

