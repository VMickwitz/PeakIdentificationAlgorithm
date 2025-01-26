function [compList, peakList, compNams, h_outs, AltMasses, AltComps, BL, dmcal, Rfit, dfvals] = fitAndAssign(fitStruct,ind,M,compPar)
%FIT_AND_ASSIGN Summary of this function goes here
%   Detailed explanation goes here

if nargin == 3
    compPar.a = 0;
end
verbate = false;
doIso = false;
dummy = false;
map = false;
unlock = false;
unlocktol = 0.10; % 0 no reassignment, 1 a lot of reassignment
dfvals = [];
nOld = 0;
param = fitStruct.param;
mz = fitStruct.mz(:);
specs = fitStruct.specs;
if isfield(fitStruct,'isoSpecs')
    doIso = true;
    isoSpecs = fitStruct.isoSpecs;
end

if isfield(param,'dummy')
    dummy = param.dummy;
    nDummy = param.n;
end

if isfield(param,'unlock')
    unlock = param.unlock;
    if isfield(param,'unlocktol')
        unlocktol = param.unlocktol;
    end
end


doOpt = false;
doAssign = true;
n = param.peakRange;
compList = [];
Hset = [];

compNams = {};
setPeaks = length(param.peakList);
AltMasses = [];
AltComps = {};

N = 0;
reassignnumber = 0;
rmpeak = 0;

A_orig = param.A;

assignDiff = M./param.R(M).*0.2;
if isfield(compPar,'assignDiff')
    assignDiff = compPar.assignDiff;
else
    compPar.assignDiff = assignDiff;
end
if isfield(compPar,'verbate')
    verbate = compPar.verbate;
end

iter = 0;
while doAssign
    iter = iter+1;
    if dummy
        if param.n == 0
            param.peakList = 0;
            break;
        end
        nDummy = param.n - length(param.peakList);
        if N == 0 && isfield(fitStruct,'peaks')
            n = param.n;
            peaks = fitStruct.peaks(:,:,param.massRange == M);
            compList = peaks(sum(~isnan(peaks),2,'omitnan') == param.n,:);
            compList(isnan(compList)) = [];
            if isempty(compList)
                peaks = peaks(sum(~isnan(peaks),2,'omitnan') > param.n,:);
                if ~isempty(peaks)
                    compList = peaks(1,1:param.n);
                else
                    warning("No suitable initial conditions for dummy fit")
                    compList = linspace(M-0.01,M+0.01,param.n);
                end
            end
        elseif N == 0
            warning("Bad initializing. Preliminary fits for fit_and_assign.")
            compList = linspace(M-0.01,M+0.01,param.n);
        elseif N == 1
            if isempty(compList)
                break;
            end
            res = sum(specs,2) - sum(ytots,2);
            [~,imax] = max(res);
            compList(end) = mz(imax);
        end

        [ytots,~,~,h_outs, par, BL] = fitPeaks(mz, specs, compList, param);
        ind = 1;
        compList = par(2:end);
        
        h_outs = h_outs(1:nDummy,:);
        if ~(nDummy == length(compList))
            error("Dummy dumb")
        end
        
        %compList
        %param.peakList
        Hcomp = sum(h_outs,2,'omitnan');
        %[Hcomp,sorter] = sort(sum(h_outs,2,'omitnan'),'descend','MissingPlacement','last');
%         sorter
%         h_outs
%         sum(h_outs,2)
%         compList
        N = 1;
    elseif N == 0
        if ind > max(n)
            param.peakList = 0;
            break;
        end
        
        [n,par,~,BL] = getFitByIndex(fitStruct,ind,M);
        
        ind = 1;
        compList = par;
        param.initBL = BL;
        [~, ~, ~, h_outs, par, BL, dmcal] = fitPeaks(mz, specs, compList, param);
        param.init_dmcal = dmcal;
        if verbate
            n
            par
            h_outs
        end

        param.initR = par(1);
        compList = par(2:end);
        Hcomp = sum(h_outs,2,'omitnan');%./param.W;
        %[Hcomp,sorter] = sort(sum(h_outs,2,'omitnan')./param.W,'descend','MissingPlacement','last');
        N = 1;
        if isnan(n)
            param.peakList = 0;
            break;
        end
        
        % Only one peak more than fitPrelim. Isotopes may reduce number of
        % peaks but increasing the number should be very rare. Also
        % significantly increases speed of assignment.
        if n(1) < param.peakRange(end)
            maxPeaks = n(1)+1;
        else
            maxPeaks = n(1);
        end
        param.peakRange = 0:maxPeaks;
        
        % If isotope signal more than 5% of total signal, refit.
        if trapz(mz,isoSpecs) > 0.05*trapz(mz,specs)
            % Continuing to next iteration of loop will cause a refit.
            compList(isnan(compList)) = [];
            iter = iter-1;
            continue;
        end
    else
        %fprintf("Mass %i refit\n",M);
        ind = 1;
        if map
            % guess_peak_number interprets A differently from other functions:
            param.A = mappingFunction(A_orig,param.W,true);
        end

        param.findA = false;
        
        if isfield(param,'tol')
            warning("Removing tol, for more flexible refitting.")
        end
        
        if unlock
            % If unlocking locked peaks allowed, then check significance of
            % locked peaks.
            dfold = dfvals;
            dfvals = evalFits(mz,specs,param);
            delta = nan(length(dfvals)-1,1);
            for i = 1:(length(dfvals)-1)
                % Difference between highest significance in newly assigned
                % peaks and significance of currently evaluated peak.
                delta(i) = max(dfvals(nOld+1:end))-dfvals(i);
            end
            % If earlier assigned peak has lost significance
            % and relative difference to old significance higher than
            % unlocktol: remove from list.
            if ~isempty(dfvals)
                indrm = (delta > dfvals(1:end-1)) & (dfvals(1:nOld) < unlocktol*dfold);
            else
                indrm = [];
            end

            if any(indrm)
                warning("Reassigned locked peak at M = %i",M)
                indrm = find(indrm,1,'first');
                reassignnumber = reassignnumber +1;
                if reassignnumber > 10
                    warning("High number of reassigned peaks. Preventing further reassignment.\n")
                else
                    peaks = param.peakList(indrm);
                    if verbate
                        fprintf("Locked peaks unlocked:")
                        peaks %#ok<NOPRT>
                    end
                    compList = [peaks, compList];
                    param.peakList(indrm) = [];
                    setPeaks = length(param.peakList);
                    compNams(indrm) = [];
                    AltComps(indrm) = [];
                    AltMasses(indrm) = [];
                    dfvals(indrm) = [];
                    
                    % Make sure to sort compList after adding the peak
                    % back, to properly initialize next guessPeakNumber
                    % call.
                    %Hcomp = [Hset(indrm); Hcomp];
                    %S = overlapScore(mz,param.peakList,compList,param.R);
                    %[~,sorter] = sort(Hcomp.*S,'descend','MissingPlacement','last');
                    [compList] = sortCompList(mz,specs,param,compList);
                    %Hcomp = Hcomp(sorter);
                    %S = S(sorter);
                    %compList = compList(sorter);
                end
            end
        nOld = setPeaks;
        end
        [n, ~, par, ~, h_outs, BL, dmcal, fitR] = guessPeakNumber(mz,specs,param,compList);

        if verbate
            compList %#ok<NOPRT>
            n %#ok<NOPRT>
            par %#ok<NOPRT>
            h_outs %#ok<NOPRT>
        end

        if isempty(n)
            peakList = [];
            h_outs = [];
            return
        end
        i1 = 1;
        if setPeaks == 0
            i1 = 0;
        end
        if any(n==0)
            i1 = i1+1;
        end
        
        
        %par(n(ind)-setPeaks+i1,:);
        compList = par(n(ind)-setPeaks+i1, setPeaks+1:end);
        param.initR = fitR(ind);
        BL = BL(ind);
        dmcal = dmcal(ind);
        Hset = sum(h_outs(1:setPeaks,:,n(ind)-min(n)+1),2);
        param.init_dmcal = dmcal;
        param.initBL = BL;

        %[Hcomp,sorter] = sort(sum(h_outs(setPeaks+1:end,:,n(ind)-min(n)+1),2,'omitnan'),'descend','MissingPlacement','last');
        % Remove sorting:

        Hcomp = sum(h_outs(setPeaks+1:end,:,n(ind)-min(n)+1),2,'omitnan');
        %sorter = 1:length(Hcomp);

        % Highest peak first is a bad approach (but better than no sorting). Reduce impact by overlap
        % measure. (H*(1-sqrt(npdf(x,peak1,width)*npdf(x,peak2,width)))) Gaussian
        % approxiamtion ok in this case
    end
    Rfit = param.initR;
    % Remove nans
    inan = isnan(compList);
    compList(inan) = [];
    Hcomp(inan) = [];
    
    if unlock
        [compList, sorter, S] = sortCompList(mz,specs,param,compList);
        Hcomp = Hcomp(sorter);
    else
        S = overlapScore(mz,param.peakList,compList,param.R);
        [~,sorter] = sort(Hcomp.*S,'descend','MissingPlacement','last');

        %Sort compList
        compList = compList(sorter);
        S = S(sorter);
        Hcomp = Hcomp(sorter);
    end

    
    %compList

    if verbate
        % Make plots if requested.
        nc = length(compList);
        np = length(param.peakList);
        Names = [string(compNams);repmat("Unassigned",[nc,1])];
        Area = [Hset;Hcomp];
        Mass = [param.peakList';compList'];
        if unlock
            dFvals = [dfvals;S];
            iterTab = table(Names,Mass,Area,dFvals);
        else
            Overlap = [nan(np,1);S];
            Score = Area.*Overlap;
            iterTab = table(Names,Mass,Area,Overlap,Score);
        end
        fprintf("Table for assignment number: %i\n\n",iter)
        disp(iterTab)
        y1s = mkFigYs(mz, [param.peakList,compList], [Hset; Hcomp], param);
        initPlot();
        plot(mz,sum(specs,2)+sum(isoSpecs,2),'k-')
        plot(mz-dmcal,sum(y1s,2)+BL(ind)+sum(isoSpecs,2),'rx')
        plot(mz,sum(isoSpecs,2),'k:')
        plot(mz-dmcal,y1s)
        for i = 1:length(param.peakList)
            xline(param.peakList(i),'k-','LineWidth',1.5)
        end
        for j = 1:length(compList)
            if j == 1
                xline(compList(j),'r--','LineWidth',1.5)
            else
                xline(compList(j),'k:','LineWidth',1.5)
            end
        end
        FWHM = mean(Mass)./param.R(mean(Mass));
        xls = [min(Mass)-2*FWHM, max(Mass)+2*FWHM];
        xlim(xls)
        title(sprintf("Assignment number: %i",iter))
        %fprintf("highest peak: %7.3f",compList(1))
    end

    if length(compList)+setPeaks < n(ind)
        compList
        setPeaks
        i1
        sorter
        par
        n
        error("CompList indexing error.")
    end
    
    if any(Hcomp == 0)
        isZero = (Hcomp == 0);
        compList(isZero) = [];
        Hcomp(isZero) = [];
        n(ind) = n(ind)-sum(isZero);
        if verbate
            y1s(:,isZero) = [];
        end
    end

    % if verbate
    %     n(ind)+setPeaks
    %     Pplot = [param.peakList,compList]
    %     Hplot = [Hset; Hcomp]
    %     %ys = y1s(:,1:length(Pplot))./Hplot(1:length(Pplot))';
    %     S = overlapScore(mz,param.peakList,compList,param.R)
    % end

    

    
    % Remove comps outside lims:
    if isfield(param,'lims')
        defect = compList-M;
        isOut = (defect < param.lims(1)) | (defect > param.lims(2));
        compList(isOut) = [];
        Hcomp(isOut) = [];
        n(ind) = n(ind) - sum(isOut);
    end
    
    %Information for peaklist_from_spectra
    if ind > max(param.peakRange)
        param.peakList = 0;
        break
    end
    
    if setPeaks == n(ind) || isempty(n)
        [isrm,indrm,dfvals,rmpeak,reassignnumber] = checkPeaks(mz,specs,param,rmpeak,reassignnumber);
        if any(indrm)
            indrm = indrm(1);
            param.peakList(indrm) = [];
            setPeaks = length(param.peakList);
            compNams(indrm) = [];
            AltComps(indrm) = [];
            AltMasses(indrm) = [];
            h_outs(indrm,:,:) = [];
            unlock = false; % Cannot allow unlocking peaks after this since nOld and fold may be mixed up.
            if isrm
                continue;
            else
                break;
            end
        else
            break;
        end
    end
    
    if ~(length(compList) == n(ind)-setPeaks) && ~dummy
        compList
        Hcomp
        length(compList)
        n
        ind
        setPeaks
        error('Faulty indexing of par')
    end
    
    looking = true;
    jnd = 1;
    % Make a copy to preserve order of compList.
    assignList = compList;
    while looking
        %[masses,compNam,doAssign,AltMass,AltComp] = assignPeaks(mz,specs,compList,param,compPar);
        [masses,compNam,doAssign,AltMass,AltComp] = assignPeaks(mz,specs,assignList,param,compPar);
        if isempty(AltMass)
            AltMass = nan;
            AltComp = "none";
        end
        
        if isempty(masses)
            %         masses
            %         compList(1)
            %         param.peakList
            %         fprintf("Breaking")

            %AltComps = [AltComps ; AltComp];
            %AltMasses = [AltMasses ; AltMass];

            if length(compList) > jnd
                jnd = jnd+1;
                %compList = [compList(2:end) compList(1)];
                assignList = [assignList(2:end) assignList(1)];
                continue;
            else
                break
            end
        elseif any(masses(1) == param.peakList)
            break
        elseif isempty(compNams)
            compNams = {char(compNam)};
            AltComps = [AltComps ; {AltComp}];
            AltMasses = [AltMasses ; {AltMass}];
        else
            compNams = [compNams ; {char(compNam)}]; %#ok<*AGROW>
            AltComps = [AltComps ; {AltComp}];
            AltMasses = [AltMasses; {AltMass}];
        end

        param.peakList = [param.peakList masses];
        setPeaks = length(param.peakList);

        if compList(jnd) ~= assignList(1)
            compList
            assignList
            jnd
            error("Check compList and assignList synchronisation.")
        end

        compList(jnd) = [];
        Hcomp(jnd) = [];
        looking = false;
        % if jnd == 1
        %     looking = false;
        % end
    end

    % if looking
    %     [isrm,indrm,dfvals,rmpeak,reassignnumber] = checkPeaks(mz,specs,param,rmpeak,reassignnumber,verbate);
    %     if any(indrm)
    %         indrm = indrm(1);
    %         param.peakList(indrm) = [];
    %         setPeaks = length(param.peakList);
    %         compNams(indrm) = [];
    %         AltComps(indrm) = [];
    %         AltMasses(indrm) = [];
    %         h_outs(indrm,:,:) = [];
    %         if isrm
    %             continue;
    %         else
    %             break;
    %         end
    %     else
    %         break;
    %     end
    % end

    if looking
        break;
    end
%     for i = 1:length(masses)
%         [~, rC] = min(abs(compList-masses(i)));
%         compList(rC) = [];
%     end
    
end

if verbate
    fprintf("Done refitting.")
end

% if dummy && ~isempty(compList)
%     comp_i = [];
%     for i = 1:length(compList)
%         res = sum(specs,2) - sum(ytots,2);
%         [~,imax] = max(res);
%         comp_i = [comp_i mz(imax)];
%         [ytots,~,~,h_outs, par, BL] = tof_fit_peaks_to_many(mz, specs, comp_i, param);
%         comp_i = par(2:end);
%         if ~(length(comp_i) == i)
%             error("Error in reinitializing final dummy fits.")
%         end
%     end
% end

j = 1;
% if ~(param.peakList == 0)
if ~(all(param.peakList == 0) && dummy)
    
    while length(compList) >= j
        
        iPL = ~(1:length(compList) == j);
        par = param;
        par.peakList = [param.peakList compList(iPL)];

        [masses,compNam,~,AltMass,AltComp] = assignPeaks(mz,specs,compList(j),param,compPar);
        
        if isempty(masses)
            if any(abs(compList(j)-param.peakList) < assignDiff/2)
                fprintf("Removing %7.4f.\n",compList(j));
                compList(j) = [];
                if verbate
                    fprintf("No comp found and within assignDiff of other peak.\n")
                end
            end
            j = j+1;
            continue;
        elseif any(masses(1) == param.peakList)
            fprintf("Removing %7.4f.\n",compList(j));
            warning("Check assignDiff. This should not be happening.")
            assignDiff
            compList(j) = [];
            continue;
        elseif isempty(compNams)
            compNams = char(compNam);
            AltComps = [AltComps ; {AltComp}];
            AltMasses = [AltMasses ; {AltMass}];
            compList(j) = [];
        else
            compNams = [compNams ; {char(compNam)}];
            AltComps = [AltComps ; {AltComp}];
            AltMasses = [AltMasses; {AltMass}];
            compList(j) = [];
        end
        fprintf("Adding mass %7.4f\n",masses(1))
        param.peakList = [param.peakList masses];
        setPeaks = length(param.peakList);
        
    end
end

% Finally remove any peaks that have zero or lower score (for some weird reason):

while true
    if ~(isempty(param.peakList) && isempty(compList))
        np = length(param.peakList);
        dfvals = evalFits(mz,specs,param,compList);
        [dfmin, indrm] = min(dfvals);

        %indrm = find(dfvals<=0);
        if dfmin <= 0%~isempty(indrm)
            if indrm > np
                compList(indrm-np)
                compList(indrm-np) = [];
            else
                irm = indrm;
                param.peakList(irm)
                param.peakList(irm) = [];
                compNams(irm) = [];
                AltComps(irm) = [];
                AltMasses(irm) = [];
                h_outs(irm,:,:) = [];
                dfvals(irm) = [];
            end
        else
            break
        end
    else
        break
    end
end

peakList = param.peakList;

if peakList == 0
    h_outs = [];
    compList = [];
    peakList = [];
    BL = 0;
else
    h_outs = squeeze(h_outs(:,:,n(ind)-min(n)+1));
    h_outs(all(isnan(h_outs(:,:)),2),:) = [];
    [~,sorter] = sort(sum(h_outs(setPeaks+1:end,:),2,'omitnan'),'descend');
    if setPeaks == 0
        % compList = compList(sorter);
        h_outs = h_outs(sorter,:);
    else
        % compList = compList(sorter);
        htemp = h_outs(setPeaks+1:end,:);
        h_outs(setPeaks+1:end,:) = htemp(sorter,:);
    end
end

if verbate
    % Make plots if requested.
    param1 = param;
    param1.peakList = [param.peakList compList];
    [~, ~, ~, h_outs, ~, BL, dmcal] = fitPeaks(mz, specs, [], param1);
    y1s = mkFigYs(mz, [param.peakList,compList], h_outs, param);
    initPlot();
    plot(mz,sum(specs,2)+sum(isoSpecs,2),'k-')
    plot(mz-dmcal,sum(y1s,2)+BL(ind)+sum(isoSpecs,2),'rx')
    plot(mz,sum(isoSpecs,2),'k:')
    plot(mz-dmcal,y1s)
    for i = 1:length(param.peakList)
        xline(param.peakList(i),'k-','LineWidth',1.5)
    end
    for j = 1:length(compList)
        xline(compList(j),'k:','LineWidth',1.5)
    end
    FWHM = mean(Mass)./param.R(mean(Mass));
    xls = [min(Mass)-2*FWHM, max(Mass)+2*FWHM];
    xlim(xls)
    title(sprintf("Final fit"))
    %fprintf("highest peak: %7.3f",compList(1))
end

end

function S = overlapScore(x,peaks,comps,Rfun)
% Returns a numeric value which roughly indicates how much comps overlaps
% with peaks. 1 = no overlap. 0 = perfect overlap.
nc = length(comps);
S = ones(nc,1);
if isempty(peaks)
    return;
end
FWHM = peaks(1)./Rfun(peaks(1));
sigma = 0.4246609*FWHM;
yp = normpdf(x,peaks,sigma);
for i = 1:nc
    [~,iclose] = min(abs(peaks-comps(i)));
    yi = normpdf(x,comps(i),sigma);
    S(i) = sqrt(1-trapz(x,sqrt(yp(:,iclose).*yi))^2);
end
end

function [isrm,indrm,dfvals,peak,rnm] = checkPeaks(mz,specs,param,rmpeak,rnm)
% Before breaking the assignment loop, check that no previously
% assigned peaks have become insignificant.
isrm = false;
indrm = false;
dfvals = [];
peak = rmpeak;
if ~isempty(param.peakList)
    % Score peak score calculated as:
    % score = numbers.*log(s*l)+A*W*s*l*fvals(1,:)
    % so if dscore = 0, and dnumbers = 1:
    % dfvals = log(s*l)/(A*W*s*l)
    n = numel(specs); % s*l
    dfmin = log(n)/(param.A*param.W*n);
    dfvals = evalFits(mz,specs,param);
    indrm = dfvals<=dfmin;
    if any(indrm)
        [~,isort] = sort(dfvals);
        indrm = indrm(isort);
        fprintf("Removing peak after assignment.\n")
        peak = param.peakList(indrm);
        isrm = true;
        rnm = rnm+1;
        if rnm > 8% || all(peak == rmpeak)
            warning("Assignment loop likely stuck in infinite loop.")
            fprintf("\n Removing peak without reassignment.\n")
            dfvals(indrm) = [];
            isrm = false;
        end
    end
end
end

function dfvals = evalFits(mz,specs,param,compList)
if nargin < 4
    compList = [];
end
testpeaks = [param.peakList compList];
testpar = param;
testpar.peakList = testpeaks;
npeaks = length(testpeaks);
dfvals = nan(npeaks,1);
[~,Fv,~,~, ~, ~] = fitPeaks(mz, specs, [], testpar);
freference = Fv(1);
ind = 1:npeaks;
for i = 1:npeaks
    testpar.peakList = testpeaks(ind~=i);
    [~,Fv] = fitPeaks(mz, specs, [], testpar);
    dfvals(i) = Fv(1)-freference;
end
end

function [compList, sorter, score] = sortCompList(mz,specs,param,compList)
nlist = length(compList);
if nlist == 0
    sorter = [];
    score = [];
    return
end
% param
% param.peakList
% compList
dfvals = evalFits(mz,specs,param,compList);
dfvals = dfvals((end-nlist+1):end);
[score, sorter] = sort(dfvals,1,'descend');
compList = compList(sorter);
end

