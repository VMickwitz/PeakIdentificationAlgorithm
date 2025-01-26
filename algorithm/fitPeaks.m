function [ytots, fvals, y1s, h_out, par, bl, dmcal] = fitPeaks(mz, specs, compList, param)
%Fit peaks of predefined peakshape and/or resolution to several spectra.
%
% In:
%     mz       - mass axis
%     specs    - matrix where each column is one of the spectrums fitted to
%     compList - initial mass list (guesses)
%     param    - structure of parameters, the fields are:
%       R:         fixed resolution value or anonymous function (R(mz))
%                  if empty it will be fitted
%       peakShape: peak shape function
%       tol:       determines how much the peak positions may deviate from
%                  the given values, can be vector with different
%                  tolerances for different peaks, default is 0.25
%       peakList:  peaks defined by mass to be included, only their height
%                  is fitted
%       RGuess:    initial value for fitting R
%       h_max:     maximum height of signals in peakList. NaN for no limit.
%       method:    Determines type of fit:
%                   1       - Fit weighted with expected errors.
%                   2 or 3  - Least squares fit.
%
% Out:
%     ytot  - matrix with the complete fitted profile for each spectrum on
%             each column (same order as input specs)
%     fvals - the first value is the value minimized by the optimisation and
%             the rest are similar numbers reflecting the accuracy of the fit to 
%             each individual spectrum.
%     y1s   - 3D matrix of y-data of all the peaks fitted to all spectra
%               1st dimension - y-data
%               2nd dimension - peak number
%               3rd dimension - spectrum number
%     h_out - areas of the fitted peaks (in mV * Th)
%     par   - the final parameters used for the fit, [R mz1 mz2 ...]
%
% 
% Sometimes fitted peak is clearly wrong (fit fails), this is indicated by the fact
% that the first value in fvals is several orders of magnitude higher than
% expected. Changing the values in compList should resolve this.

mz = mz(:);
specs = specs(:,:);
compList = compList(:)';
compList(isnan(compList)) = [];
param.fitR = false; % fit resolution
Rtol = false;
tol = 0.1;
useTol = false;
param.usePeakList = false;
initR = 5000;
param.fitBaseline = true;
param.simple = true; % fits one baseline to all factors
param.deltaMCal = false;
param.dmc = 0;

% interpret arguments

if nargin == 4
    if isfield(param, 'peakList')
        param.peakList = param.peakList(:)';
        param.peakList(isnan(param.peakList)) = [];
        if ~isempty(param.peakList)
            param.usePeakList = true;
        end
    end
    if isfield(param, 'R')
        if isempty(param.R)
            param.fitR = true;
            param.Rtol = [];
        elseif isa(param.R, 'function_handle')
            if ~isempty(compList)
                param.R = param.R(mean(compList));
            elseif param.usePeakList
                param.R = param.R(mean(param.peakList));
            else
                param.R = param.R(mean(mz));
            end
        end
        if isfield(param,'Rtol')
            % tolerance for variation in width (percent)
            if param.usePeakList && param.Rtol > 0
                param.fitR = true;
                Rtol = true;
                Rlim = 1+[1 0 -1].*param.Rtol.*1e-2;
                Rlim = mean(mz)./((mean(mz)./param.R).*Rlim);
            end
        end
    else
        param.fitR = true;
    end
    if isfield(param, 'RGuess') && param.fitR
        initR = param.RGuess;
    end
    if isfield(param, 'tol')
        tol = param.tol(:)';
        useTol = true;
    end
    if isfield(param, 'peakShape') && ~isfield(param.peakShape, 'sp')
        pp = spline(param.peakShape.dat(:,1),param.peakShape.dat(:,2));
        param.peakShape.sp = pp;
    elseif ~isfield(param, 'peakShape')
        rxi = -30:0.1:30;
        param.peakShape.sp = spline(rxi, normpdf(rxi, 0, 1));
    end
    if ~isfield(param,'A')
        param.A = 1;
    end
    if ~isfield(param,'weights')
       param.weights = ones(1,size(specs,2));
    end
    if isfield(param,'BL') && ~isempty(param.BL)
        if length(param.BL) > 1
            error("Baseline input too long. Select correct value before fitting.")
        end
        bl = param.BL;
        param.fitBaseline = false;
        if isfield(param,'baseline')
            param.baseline = param.baseline+bl;
        else
            param.baseline = ones(size(specs))*bl;
        end
    else
        bl = 0;
    end
    if ~isfield(param,'baseline')
        param.baseline = 0;
    end
    if isfield(param,'dmcal') && param.dmcal > 0
        % tolerance for variation in mass calibration (ppm)
        if param.usePeakList
            if isfield(param,'dmfit')
                param.dmc = param.dmfit;
            else
                param.deltaMCal = true;
                dmcal = mean(mz).*([-1 1].*param.dmcal*1e-6);
                if isfield(param,'init_dmcal')
                    init_dmcal = param.init_dmcal;
                else
                    init_dmcal = 0;
                end
            end
        end
    end
    if ~isfield(param,'method')
        % 1 = chi^2 fit
        % 2 = MSE
        param.method = 1;
    end
else
    error("Parameters not an optional arguemnt.")
end


% set values and fit
init = [];
ll = [];
ul = [];
param.calc = false;
nrP = length(compList);

if length(tol) > nrP
    tol = tol(1:nrP);
end

n = size(specs,2);

warning('error','MATLAB:rankDeficientMatrix');

if nrP == 0 && ~param.usePeakList
    param.fitR = false;
end

% Set initial values and limits:
if nrP > 0
    init = compList;
    if useTol
        % Prevent overlap:
        ioverlap = false;
        if ~isempty(compList)
            [sortCL, isort] = sort(compList);
            dCL = diff(sortCL);
            ioverlap = diff(sortCL) < 2*tol;
            limOL = sortCL(1:end-1) + dCL/2;
        end
        ll = compList - tol;
        ul = compList + tol;
        if any(ioverlap)
            iOL = find(ioverlap);
            indLL = isort(iOL+1);
            indUL = isort(iOL);
            ll(indLL) = limOL(iOL);
            ul(indUL) = limOL(iOL);
        end
        % Prevent crossing over locked peak locations:
        for i = 1:length(param.peakList)
            pli = param.peakList(i);
            iCorr = ll < pli & ul > pli;
            if any(iCorr)
                ill = compList > pli & iCorr;
                iul = ~ill & iCorr;
                ll(ill) = param.peakList(i);
                ul(iul) = param.peakList(i);
            end
        end
    else
        % If no tol specified, allows peak anywhere within mass axis limits
        ll = ones(1,nrP).*mz(1);
        ul = ones(1,nrP).*mz(end);
    end
    if any(init <= ll)
        ind = init<=ll;
        init(ind) = ll(ind)+(ul(ind)-ll(ind))/10;
    end
    if any(init >= ul)
        ind = init>=ul;
        init(ind) = ul(ind)-(ul(ind)-ll(ind))/10;
    end
end

if param.deltaMCal
    init = [init_dmcal init];
    ll = [dmcal(1) ll];
    ul = [dmcal(2) ul];
end

if param.fitBaseline
    if isfield(param,'initBL')
        initBL = double(param.initBL);
    else
        initBL = mean(0.5*(mean(specs(1:3,:)+specs(end-2:end,:)))./param.weights);
    end
    if param.simple
        try
            init = [initBL init];
        catch ME
            init
            initBL
            rethrow ME
        end
        if initBL > 0
            ll = [0.5*initBL ll];
            ul = [2*initBL ul];
        else
            warning("Fitting negative baseline.")
            ll = [2*initBL ll];
            ul = [0 ul];
        end
    else
        init = [(specs(1,:)+specs(end,:))/2 init];
        ll = [zeros(1,n) ll];
        ul = [specs(1,:)+specs(end,:) ul];
    end
end

if param.fitR
    if Rtol
        init = [Rlim(2) init];
        ll = [Rlim(1) ll];
        ul = [Rlim(3) ul];
    end
    if isfield(param,'initR')
        init(1) = param.initR;
    end
end

init = double(init);
ll = double(ll);
ul = double(ul);

% Make the fit
if nrP > 0 || param.fitR || param.fitBaseline
    opt = optimset('display', 'off', 'TolX', 1e-5, 'TolFun', 1e-16,'MaxFunEvals',100*(nrP+1));
    par = boundfmin(@(par) fitfun(par, mz, specs-bl, param), init, ll, ul, opt);
    if ~param.fitR
        par = [param.R par];
    end
else
    par = [param.R bl];
    param.fitBaseline = true;
end

% calculate the final parameters
param.calc = true;
param.fitR = true; % par has already been converted to the form used when fitting resolution

[fval, fvals, ytots, h_out, kernel] = fitfun(par, mz, specs, param); % use fitfun to calculate the final parameters

%h_out = h_out(:, 1);

if param.fitBaseline
    if param.simple
        bl = par(2);
        par(2) = [];
    else
        bl = par(2:n+1);
        par(2:n+1) = [];
    end
    ytots = ytots + bl.*param.weights;
end

if param.deltaMCal
    dmcal = par(2);
    par(2) = [];
else
    dmcal = param.dmc;
end

if param.usePeakList
    nrP = nrP + length(param.peakList);
end

l = size(kernel,1);
y1s = NaN(l,nrP,n);

for j = 1:n
    y1 = NaN(l,nrP);
    for i = 1:nrP
        y1(:, i) = kernel(:, i) * h_out(i,j);
    end
    y1s(:,:,j) = y1;
end

fvals = [fval; fvals]; %first element in fvals represents total error as determined by fitfun.

end
%%
function [fval, fvals, yhats, h_outs, kernel] = fitfun(par, x, ys, param)

if param.fitR
    param.R = par(1);
    par(1) = [];
end
if param.fitBaseline
    if param.simple
        ys = ys-par(1).*param.weights;
        param.baseline = param.baseline + par(1).*param.weights;
        par(1) = [];
    else
        n = size(ys,2);
        ys = ys-par(1:n);
        param.baseline = param.baseline + par(1:n);
        par(1:n) = [];
    end
end
if param.usePeakList
    par = [par param.peakList];
end
if param.deltaMCal
    dmcal = par(1);
    par(1) = [];
    par = par-dmcal;
else
    par = par-param.dmc;
end

[fval, fvals, yhats, h_outs, kernel] = calcFit(par, x, ys, param);
    
end
%%
function [fval, fvals, yhats, h_outs, kernel] = calcFit(par, x, ys, param)%, h_lim)

% This function is similar to local function fitfun in tof_fit_peaks.m but
% works for several functions. Large portions of the code is copied from
% this function. Horribly poorly optimized. Improve for performance.

n = size(ys,2); %number of spectra fitted to
l = size(ys,1); %length of data fitted to

nrP = length(par);

kernel = zeros(length(x), nrP);

for i = 1:nrP
    % fast kernel generation without checks
    
    w = 0.4246609 * round(par(i)) / param.R;
    xi = (x - par(i)) / w;

    index = discretize(xi, [-inf, param.peakShape.sp.breaks(2:param.peakShape.sp.pieces), inf]);

    xs = xi - param.peakShape.sp.breaks(index)';

    v = param.peakShape.sp.coefs(index, 1);
    for j = 2:param.peakShape.sp.order
        v = xs(:) .* v + param.peakShape.sp.coefs(index, j);
    end
    
    kernel(:, i) = max(v, 0);
    
    if param.calc
        if ~all(kernel(:,i)==0)
            kernel(:, i) = kernel(:, i) / trapz(x, kernel(:, i));
        end
    end
end

if any(sum(kernel == kernel(1,:), 1) >= (l-1))
    % If attempting to fit peaks too far to one side, there may be constant
    % kernels -> rank deficiency problems.
    %irm = all(kernel == kernel(1,:), 1);
    irm = sum(kernel == kernel(1,:), 1) >= (l-1);
    kernel(:,irm) = [];
else
    irm = 0;
end

if nrP == 0
    h_outs = [];
else
    try
        h_outs = kernel \ ys;
    catch
        [~,ih,ik] = unique(par(~irm));

        kernel2 = kernel(:,ih);
        %[kernel2,ih,ik] = unique(kernel','rows');

        if any(all(round(kernel2 - kernel2(1,:),5)==0, 1))
            % If attempting to fit peaks too far to one side, there may be
            % constant kernels -> rank deficiency problems. Remove them
            i0 = all(round(kernel2 - kernel2(1,:),5)==0, 1);
            kernel2(:,i0) = [];
        else
            i0 = 0;
        end

        try
            h_outs = kernel2 \ ys;
        catch
            [~,ih,ik] = unique(round(par(~irm),3));
            kernel2 = kernel(:,ih);

            if any(all(round(kernel2 - kernel2(1,:),4)==0, 1))
                % If attempting to fit peaks too far to one side, there may be
                % constant kernels -> rank deficiency problems. Remove them
                i0 = all(round(kernel2 - kernel2(1,:),4)==0, 1);
                kernel2(:,i0) = [];
            else
                i0 = 0;
            end

            try
                h_outs = kernel2 \ ys;
            catch ME
                warning("Rank deficiency problem (several linearly dependent kernels).")
                warning("Attempting to ignore")
                warning(ME.message)
                h_outs = zeros(size(kernel2,2),size(ys,2));
            end
        end

        if any(i0)
            % Add i0 values back to h_outs
            h_temp = zeros(size(kernel2,2)+sum(i0),size(ys,2));
            h_temp(~i0,:) = h_outs;
            h_outs = h_temp;
        end

        h_outs = h_outs(ik,:);
        ih = setdiff(1:(nrP-sum(irm)),ih);
        h_outs(ih,:) = 0;
    end
end

if any(irm)
    % Add irm values back to h_outs and kernel
    h_temp = zeros(nrP,size(ys,2));
    h_temp(~irm,:) = h_outs;
    h_outs = h_temp;

    k_temp = zeros(size(ys,1),nrP);
    k_temp(:,~irm) = kernel;
    kernel = k_temp;
end

if any(any(h_outs(:, :) < 0))
    t.Display = 'off';
    for i = 1:n
        if any(h_outs(:,i) < 0)
            try
                [h_out, ~, ~, ~, ~] = lsqnonneg(kernel, double(ys(:,i)), t);
                h_outs(:,i) = h_out(:,1); %#ok<AGROW>
            catch
                % Catch rank deficiency problems (fix this)
                [~,ih,ik] = unique(par(~irm));
                kernel2 = kernel(:,ih);
                [h_out, ~, ~, ~, ~] = lsqnonneg(kernel2, double(ys(:,i)), t);
                h_outs(:,i) = h_out(ik,1); %#ok<AGROW>
                ih = setdiff(1:nrP,ih);
                h_outs(ih,i) = 0; %#ok<AGROW>
            end
        end
    end
end

if nrP == 0
    yhats = zeros(size(ys));
else
    yhats = kernel * h_outs;
end

if isfield(param,'baseline') && param.method == 1
    if any(yhats+param.baseline == 0)
        param.baseline(yhats+param.baseline == 0) = realmin;
    end
    fvals = sum(((ys-yhats).^2)./(yhats+param.baseline),1,'omitnan')';
    fval = sum(fvals,'omitnan')./(l*n);
else
    fvals = vecnorm((ys-yhats),2,1)';
    fval = sum(fvals.^2,'omitnan');
    if param.method == 1
        warning('No baseline correction value given.')
    end
end

end

