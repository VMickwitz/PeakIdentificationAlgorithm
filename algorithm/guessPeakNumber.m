function [numbers, score, peaks, removed, h_outs, bls, dmcal, Rfits] = guessPeakNumber(mz, specs, param, initList)
% Function for determining the amount of peaks to fit to a spectrum.

% Estimates the true number of distributions a peak is composed of
% Inputs: 
%       mz - mass axis
%       spec - the spectra distributions are fitted to. Matrix where each
%       column is one spectrum
%       param - structure of parameters with the following fields
%           -R : Fixed resolution or anonymous function R(mz). Fitted if
%           empty (NOT RECOMMENDED)
%           -peakShape : peak shape function. Assumes normal distribution
%           if missing
%           -peakList : Masses of known peaks in spectrum.
%           -peakRange : Array of possible numbers of peaks (guess).
%           Assumes between 1 and 6 if not given.
%
% Outputs: 
%       numbers - array containing the probable amounts of peaks in the
%       spectra
%       score - Parameters describing roughly how likely the corresponding
%       number of peaks is.
%       peaks - positions of peaks in the likely amounts of peaks, column i
%       representing peak number i.
%
%

mz = mz(:);
specsfit = specs(:,:);
useR = true;
specTot = sum(specsfit,2);
[~, iMax] = max(specTot);
peakL = 0;
weights = 1;
useList = false;
listL = 0;

% Processing inputs:
if nargin < 2
    error('Too few input arguments')
elseif nargin == 2
    warning('Estimating peaknumbers without known resolution is not recommended')
    useR = false;
elseif nargin > 2
    if isfield(param,'peakRange')
        peakRange = param.peakRange;
        peakRange = sort(peakRange);
        if peakRange(1) == 0
            peakRange = 0:max(peakRange);
        else
            peakRange = 1:max(peakRange);
        end
    else
        warning("using default peakRange.")
        peakRange = 1:6;
    end   
    if ~isfield(param,'R')
        warning('Estimating peaknumbers without known resolution is not recommended')
        useR = false;
    elseif isa(param.R, 'function_handle')
        param.R = param.R(mz(iMax));
    end
    if isfield(param, 'peakShape') && ~isfield(param.peakShape, 'sp')
        pp = spline(param.peakShape.dat(:,1),param.peakShape.dat(:,2));
        param.peakShape.sp = pp;
    elseif ~isfield(param, 'peakShape')
        rxi = -30:0.1:30;
        param.peakShape.sp = spline(rxi, normpdf(rxi, 0, 1));
    end
    if isfield(param,'peakList')    
        if isempty(param.peakList)
            param = rmfield(param,'peakList');
        else
            peakL = length(param.peakList);
        end
    end
    if isfield(param,'weights')
        weights = param.weights;
        %specsfit = specsfit.*weights;
    end
    if ~isfield(param,'peakThreshold')
        param.peakThreshold = 0;
    end
    if ~isfield(param,'doRm')
        param.doRm = false;
    end
    if isfield(param,'W')
        Wa = param.W;
    else
        Wa = 1;
    end
end

if nargin > 3 && ~isempty(initList)
    % Init list is expected to be sorted in order of likelyhood beforehand
    useList = true;
    listL = length(initList);
end

peakRange(peakRange < peakL) = [];

numbers = peakRange;
a = size(peakRange,2);


% Check for very noisy spectra by comparing sides to height
% Spectrum modified if Height is less than 0.01 times both sides
remove = [];
rm = false;

s = size(specsfit,2);
% testArea = mean(trapz(mz,specsfit))

%Height limit for data not analysed:
if param.peakThreshold > 0
    limSmall = param.peakThreshold/s;
    
    for i=1:s
        H = max(specsfit(:,i));
        
        if H < limSmall
            rm = true;
            remove = [remove i]; %#ok<AGROW>
        elseif  H < (20*specsfit(1,i)) && H < (20*specsfit(end,i))
            %remove = [remove i]; %#ok<AGROW>
            %rm = true;
            specsfit(:,i) = specsfit(:,i)-((specsfit(1:2,i))+(specsfit(end-1:end,i)))/4;
%         elseif any(specsfit < 0,'all')
%             if any(specsfit() < -H/20,'all')
%                 rm = true;
%                 remove = [remove i];
%             end
        end
    end
    specsfit(specsfit<0) = 0;
end

%Check for skinny peaks (no peak should have higher R than expected value)
if useR && param.doRm
    for i=1:s
        param1 = rmfield(param,'R');
        if isfield(param1,'peakList')
            param1 = rmfield(param1,'peakList');
        end
        [~, p] = max(specsfit(:,i));
        [~, ~, ~, ~, par] = fitPeaks(mz, specsfit(:,i), mz(p), param1);
        if par(1) > (1.1*param.R) && ~any(remove == i)
            remove = [remove i]; %#ok<AGROW>
            rm = true;
        end
        %change this to comparing fitted R with m/FWHM?
        %--> trouble with double peaks?
    end
    
end

if rm
    if length(remove) == s
        warning('No factor deemed fit for analysis')
        numbers = [];
        score = [];
        removed = [];
        peaks = [];
        h_outs = [];
        return
    end
    specsfit(:,remove) = [];
    s = size(specsfit,2);
    %warning('Spectrum/spectra were removed from optimization')
end

removed = remove;

specTot = sum(specsfit,2);
[H, iMax] = max(specTot);
As = trapz(mz,specsfit,1);
A = sum(As);
W = A/(H); %estimate of peak width for generating compList
w = round(W / (mz(2) - mz(1))); %assumes evenly spaced x-axis

% Preparing arrays:
l = size(mz,1);
fvals = nan(s+1,a);
pars = nan(a,numbers(a)+1);
ytots = nan(l,s);
bls = zeros(a,1);
%resTot = nan(a,1);
dmcal = nan(a,1); % Currently only one dmcal fit ofr all sepctra.
%y1s = NaN(l,a,s);

compList = nan(a-peakL,1);

for i = 1:a
    %compList = genCompList(peakRange(i),mz,w,iMax);
    if useList && (listL >= (numbers(i)-peakL))
        if i == 1
            if peakL == 0
                if numbers(1) == 0
                    %compList = [];
                else
                    compList(numbers(1)) = initList(1);
                end
            elseif numbers(1) ~= peakL
                compList(1:(numbers(1)-peakL)) = initList(1:(numbers(1)-peakL));
            end
        elseif numbers(i) ~= peakL
            %compList(numbers(i)-peakL) = initList(numbers(i)-peakL);
            compList(1:(numbers(i)-peakL)) = initList(1:(numbers(i)-peakL));
        % else
        %     compList(numbers(i)-peakL) = initList(numbers(i)-peakL);
        end
    elseif i == 1 && numbers(i) - peakL == 1
        compList(1) = mz(iMax);
    elseif i==1 && numbers(i) - peakL > 1
        compList(1:numbers(1)-peakL) = genCompList(numbers(1)-peakL,mz,w,iMax);
    elseif i == 1 && numbers(1) == 0
        %compList = [];
    elseif ~(i == 1)

        % maximum of residual:
        %[~, ind] = max(sum(specsfit-ytots,2),[],1);
        if numbers(i) == 0
            %compList = [];
            error("This should never happen")
        else
            % Weighted maximum of smoothed residual:
            k = round(mz(iMax)./(param.R*(mz(2)-mz(1)))); % Smoothing window ~ FWHM
            % Force edges to zero
            smoothRes = smooth([0; 0; sum((specsfit-ytots)./sqrt(ytots+param.baseline+2*bl.*param.weights),2); 0; 0],k);
            %smoothRes = smooth([0; 0; abs(sum((specsfit-ytots)./sqrt(ytots+param.baseline+2*bl.*param.weights),2)); 0; 0],k);
            [~, ind] = max(smoothRes(4:end-3),[],1);
            ind = ind+1; % Account for edges not being included
            % initPlot();
            % plot(mz,smoothRes(3:end-2))

            %[~, sMax] = max(fvals(2:end,(i-1)));
            %[~, ind] = max(specsfit(:,sMax)-ytots(:,sMax));
            compList(numbers(i)-peakL) = mz(ind);
        end
    end
    if i > 1
        param.initBL = bls(i-1,:);
    end
    %h_outs(1:numbers(i),:,i)
    %compList
    % numbers(i)
    % round(mz(1))
    % % i
    % [ytots(:,:) , fvals(:,i), ~, ~, pars(i,1+peakL:numbers(i)+1),bl] = tof_fit_peaks_to_many(mz, specsfit, compList, param);

    % fprintf("Guess initialization.")
    % numbers(i)
    % compList
    [ytots(:,:) , fvals(:,i), ~, ~, pars(i,1+peakL:numbers(i)+1),bl,dmcal(i)] = fitPeaks(mz, specsfit, compList, param);

    % par contains residual which is removed later.
    
    if isfield(param,'peakList')
        pars(i,1:peakL) = param.peakList;
    end
    
    if ~useList
        if ~(i == 1) && (fvals(1,i) > fvals(1,i-1))    %Checking that fitting was successful
            compList2 = genCompList(numbers(i)-peakL,mz,w,iMax);
            [ytots(:,:) , fvals(:,i), ~, ~, pars(i,1+peakL:(numbers(i)+1)),bl,dmcal(i)] = fitPeaks(mz, specsfit, compList2, param);
            if fvals(1,i) > fvals(1,i-1)
                %warning('unsuccessful fit at peak number %i',numbers(i))
                doStart = false;
            else
                compList = compList2; % if list 2 better than one, use that.
            end
        end
    end

    %resTot(i) = mean(trapz(mz,abs(specsfit-ytots)));

    if ~(numbers(i) == 0)
        compList(1:numbers(i)-peakL) = pars(i,2+peakL:(numbers(i)+1));
    end
    bls(i,:) = bl;
    
end
Rfits = pars(:,1+peakL);
pars(:,1+peakL) = []; %remove R from peaks
% numbers
% pars

if isempty(fvals)
    error('Major error in fitting')
end

fvals(2:end,:) = fvals(2:end,:).^2;

if param.findA
    param1 = param;
    param1.minSignal = max(max(sum(specs,2))/100,param.minSignal);
    param.A = findA(mz,s,param1);
    fprintf("A = %4.2f \n",param.A)
    scores = numbers.*log(s*l)+param.A*s*l*fvals(1,:);
elseif isfield(param,"A")
    scores = numbers.*log(s*l)+param.A*Wa*s*l*fvals(1,:);
else
    scores = numbers.*log(s*l)+s*l*fvals(1,:);
end
% For efficiency, limit the number of returned fits.
start = 1;
stop = max(peakRange)-numbers(1)+1;


%calculate final h_outs
h_outs = nan(numbers(stop),size(specs,2),stop-start+1);
inc = 1:size(specs,2);
inc(remove) = [];

param1 = param;
if isfield(param,'Rtol')
    param1 = rmfield(param1,'Rtol');
end

for i = start:stop
    param1.peakList = pars(i,1:end);
    param1.peakList(isnan(param1.peakList)) = [];
    param1.baseline = param.baseline;
    param1.BL = bls(i,1);
    param1.R = Rfits(i,1);
    param1.dmcal = 0;
    [~,~,~,h_outs(1:length(param1.peakList),inc,i-start+1),~,~] = fitPeaks(mz, specsfit, [], param1);
end

numbers = numbers(start:stop);
scores = scores(start:stop);
bls = bls(start:stop);
dmcal = dmcal(start:stop);
Rfits = Rfits(start:stop);
[score, I] = sort(scores);
numbers = numbers(I);
bls = bls(I,1);
dmcal = dmcal(I,:);
Rfits = Rfits(I,:);
nanS = isnan(score);
score(nanS) = [];
numbers(nanS) = [];
peaks = pars;
end
%%
function [compList] = genCompList(n, x, w, iMax)
%Generates (very rough) initial values for fitting function
%   n - number of values
%   w - number of indexes in width of peak
%   iMax - index of highest point on peak

x = x(2:end-1); % avoid initializing at edges
xlen = size(x,1);
if n == 1
    indexes = iMax;
else
    w = round(w/(2*(n-1)))*2*(n-1);
    start = iMax - w/2;
    
    if start+w > xlen
        start = xlen-w;
    end
    if start < 1
        w = w + 1-start;
        start = 1;
    end
    indexes = NaN(n,1);
    for i = 1:n
        indexes(i) = start + (i-1)*w/(n-1);
    end
end
try
    compList = x(indexes);
catch
    compList = ones([1 n]).*x(iMax);
    warning("Error in genCompList>guess_peak_number")
end
end
