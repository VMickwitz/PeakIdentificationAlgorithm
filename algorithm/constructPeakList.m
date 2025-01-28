function [finalStructure] = constructPeakList(fitStruct)
%CONSTRUCTPEAKLIST Summary of this function goes here
%   Detailed explanation goes here

param = fitStruct.param;
massRange = param.massRange;
mz = fitStruct.mz;
specs = fitStruct.specs;
M = round(mz);
peakThreshold = 0;
compPar = struct;
dummy = false;
useCompF = false;
doIsotopes = true;

if isfield(param,'BL')
    param = rmfield(param,'BL');
end
if isfield(fitStruct,'compPar')
    compPar = fitStruct.compPar;
end
if isfield(param,'compF')
    useCompF = true;
    if isstruct(param.compF)
        % if compF is a list, reformat into function:
        compF = param.compF;
        param.compF = @(x,par) findComp(x,param.compF,par);
    % elseif isa(param.compF,'function_handle')
    %     param.compF = @(x,par) param.compF(x,par);
    end
else
    error("No compF provided.")
end
if ~isfield(param,'peakList')
    param.peakList = [];
end
if isfield(param,'Isotopes')
    doIsotopes = fitStruct.param.Isotopes;
end
if isfield(param,'peakThreshold')
    peakThreshold = param.peakThreshold;
end
if isfield(param,'dummy')
    dummy = param.dummy;
end

peakList = param.peakList;

NamAll = {};
MasAll = {};
AltComps = {};
AltMasses = [];
isoNam = {};
isoPeaks = {};
compNam = {};
relH = [];
significance = [];

isotopes = [];
broken = 0;
h = waitbar(0,'Making peak list...');
steps = massRange(end)-massRange(1)+1;

n = length(massRange);
pmax = max(param.peakRange);
s = size(specs,2);
fit.param = param;
fit.peaks = nan(pmax,n);
fit.H = nan(pmax,s,n);
fit.numbers = nan(n,1);
fit.names = strings(pmax,n);
BL = nan(n,1);
dmcal = nan(n,1);
Rfit = nan(n,1);

if length(param.A) == 1
    param.A = ones(n,1)*param.A;
end

fitStruct.param = param;
PL = round(peakList);
% Construct peakList based on fits:
fitIso = fitStruct;
fitIso.param.BL = zeros(n,1);
isoSpecs = zeros(size(specs));

for i = 1:n
    M_i = massRange(i);

    list.compNams = NamAll;
    list.peakList = cell2mat(MasAll);
    
    
    %prepare peakList
    param.peakList = peakList(PL==M_i);
    peaks = [];
    
    i_M = M==M_i;
    ind = 1;

    if doIsotopes
        isoY = getIsoYs(fitIso,list,M_i);
    else
        isoY = zeros(sum(i_M),s);
    end

    specs(i_M,:) = specs(i_M,:) - isoY;
    fitIso.specs = specs;
    fitStruct.param.baseline(i_M,:) = fitStruct.param.baseline(i_M,:) + isoY;
    fitIso.param.baseline(i_M,:) = fitStruct.param.baseline(i_M,:);
    isoSpecs(i_M,:) = isoSpecs(i_M,:) + isoY;
    fit_i = fitStruct;
    fit_i.param.W = param.W(i);
    fit_i.param.A = param.A(i);
    fit_i.param.baseline = fitStruct.param.baseline(i_M,:)./param.W(i);
    fit_i.mz = fitStruct.mz(i_M);
    %A1 = mean(trapz(fit_i.mz,specs(i_M,:)).*param.weights)
    fit_i.specs = specs(i_M,:)./param.W(i);
    fit_i.isoSpecs = isoY./param.W(i);
    %A2 = mean(trapz(fit_i.mz,fit_i.specs).*param.weights)

    %Skip tiny signals
    if max(fit_i.specs.*param.W(i),[],'all') < peakThreshold
        continue
    end
    
    if dummy
        fit_i.param.n = param.n(i);
    end
    
    [comps, peaks, nams, ~, AltM, AltC, BLinit, dminit, Rinit, dfvals] = fitAndAssign(fit_i,1,M_i,compPar);
    
    fit_i.param.init_dmcal = dminit;
    fit_i.param.initBL = BLinit;
    fit_i.param.initR = Rinit;
    
    if dummy && isempty(comps) && isempty(peaks)
        continue;
    end

    [h_outs,BL(i),dmcal(i),Rfit(i)] = getHfit(fit_i,[peaks comps],M_i);
    
    if size(h_outs,1) ~= (length(peaks)+length(comps))
        warning("Critical error in fitted heights. Adding zeros.")
        h_outs
        peaks
        comps
        dn = (length(peaks)+length(comps)) - size(h_outs,1);
        h_outs = [h_outs; zeros(dn,size(h_outs,2))];
    end
%     h_outs = h_outs.*permute(param.W(i),[2 3 1]);
    
%     
    h_set = sum(h_outs,2,'omitnan');
    

    if ~isempty(comps)
        for j = 1:length(comps)
            compNam = [compNam; "unknown"]; %#ok<AGROW>
            if useCompF
                [NamS,compMass] = param.compF(comps(j),compPar);
            end
            if isempty(compMass)
                AltM = [AltM; {nan}]; %#ok<AGROW>
                AltC = [AltC; {"none"}]; %#ok<AGROW>
            else
                AltM = [AltM; {compMass}]; %#ok<AGROW>
                AltC = [AltC; {NamS}]; %#ok<AGROW>
            end
        end
    else
        compNam = [];
    end
    
    if isempty(nams)
        nams = [];
    end
    % Make this work for several indices
    % if isempty(AltM)
    %     AltM = nan;
    %     AltC = {"none"};
    % end
    try
        relH = [relH; h_set./sum(h_set,1,'omitnan')];
        NamAll = [NamAll; isoNam; nams; compNam];
        MasAll = [MasAll; isoPeaks; num2cell(double(peaks')); num2cell(double(comps'))];
        AltComps = [AltComps; AltC];
        AltMasses = [AltMasses; AltM];
        significance = [significance;dfvals];

        pn = length(peaks)+length(comps);
        fit.peaks(1:pn,i) = [peaks,comps];
        fit.H(1:pn,:,i) = h_outs*param.W(i); % Scale to original spectrum
        fit.numbers(i) = pn;
        fit.names(1:pn,i) = [string(nams)', compNam'];
    catch me
        pn
        peaks
        comps
        nams
        compNam
        NamAll
        isoNam
        rethrow(me)
    end
    
    isoNam = {};
    isoPeaks = {};
    compNam = {};

    waitbar((M_i+1-massRange(1))/steps)
end
try
    close(h)
catch
    F = findall(0,'type','figure','tag','TMWWaitbar');
    delete(F);
end
specsFinal = specs;
PeakList = cell2mat(MasAll);
CompNams = string(NamAll);
AltList = AltMasses;
AltNams = AltComps;

finalStructure = fitStruct;
finalStructure.H = fit.H;
finalStructure.numbers = fit.numbers;
finalStructure.peaks = fit.peaks;
finalStructure.PeakList = PeakList;
finalStructure.CompNams = CompNams;
finalStructure.Significance = significance;
finalStructure.AltList = AltList;
finalStructure.AltNams = AltNams;
finalStructure.relH = relH;
finalStructure.specsFinal = specsFinal;
finalStructure.mz = mz;
finalStructure.param = fitStruct.param;
finalStructure.param.BL = BL.*param.W; % Scale to original spectrum
finalStructure.isoSpecs = isoSpecs;
finalStructure.deltaMassCal = dmcal;
finalStructure.Rfit = Rfit;
if useCompF
    finalStructure.param.compF = compF;
end

end

