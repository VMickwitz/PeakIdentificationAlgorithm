function [y1s] = mkFigYs(x, peaks, h_in, param)
%MK_FIG_DATA Summary of this function goes here
%   Detailed explanation goes here

% Create mass axis
if ~isfield(param,'massRange')
    param.massRange = [1 n];
end

if ~isfield(param,'R')
    if isfield(param,'fwhm')
        param.R = @(m) m/param.fwhm;
    else
        warning('No given resolution, using default value of fwhm = 0.1')
        param.R = @(m) m/0.1;
    end
end
if ~isfield(param,'W')
    param.W = 1;
end
if ~isfield(param,'H_scale')
    param.H_scale = 1;
end

if ~isfield(param,'baseline')
    param.baseline = param.noise;
end
if ~isfield(param,'Area')
    param.Area = true;
    % If true, peak size parameters correspond to area of peak, not height.
end
if isfield(param, 'peakShape') && ~isfield(param.peakShape, 'sp')
    pp = spline(param.peakShape.dat(:,1),param.peakShape.dat(:,2));
    param.peakShape.sp = pp;
elseif ~isfield(param, 'peakShape')
    rxi = -30:0.1:30;
    param.peakShape.sp = spline(rxi, normpdf(rxi, 0, 1));
end

if isrow(x)
    x = x';
end

if any(isnan(peaks))
    inan = isnan(peaks);
    peaks(inan) = [];
    h_in(inan,:) = [];
end

y1s = make_peaks(x,peaks,h_in,param.R,param.peakShape,param.Area);

function [y_peaks] = make_peaks(x, par, H, R, peakShape, isA)

nrP = length(par);
kernel = zeros(length(x), nrP);

for i = 1:nrP
    % fast kernel generation without checks

    w = 0.4246609 * round(par(i)) / R(par(i));
    xi = (x - par(i)) / w;
    
    % custom pp-evaluation
    %[~, index] = histc(xi, [-inf, peakShape.sp.breaks(2:peakShape.sp.pieces), inf]);
    
    index = discretize(xi, [-inf, peakShape.sp.breaks(2:peakShape.sp.pieces), inf]);
    xs = xi - peakShape.sp.breaks(index)';
    v = peakShape.sp.coefs(index, 1);
    for j = 2:peakShape.sp.order
        v = xs(:) .* v + peakShape.sp.coefs(index, j);
    end
    
    kernel(:, i) = max(v, 0);
    
    if isA
        kernel(:, i) = kernel(:, i) / trapz(x, kernel(:, i));
    elseif max(kernel(:, i)) > 0
        kernel(:, i) = kernel(:, i) / max(kernel(:, i));
    else
        kernel(:,i) = zeros(length(x),1);
    end
end

y_peaks = permute(kernel,[2 3 1]) .* H;
y_peaks = permute(y_peaks,[3 1 2]);
end
end

