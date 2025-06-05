function [specsAll,mzAll,param] = importDataset(dataPath,fnames,param)
%IMPORTDATASET Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2 || isempty(fnames)
    files = dir(dataPath);
    files(1:2) = [];
    fnames = cell(size(files));
    [fnames{:}] = files.name;
end
if nargin < 3
    param = struct;
end

isGUI = false;
doTime = false;
doRange = false;
isProcessed = false;

if isfield(param,'guiBar')
    isGUI = true;
    param.guiBar.Message = "Importing Spectra...";
    param.guiBar.Value = 0;
else
    h = waitbar(0,"Importing spectra...");
end
if isfield(param,'time') && param.time
    doTime = true;
end
if isfield(param,'massRange')
    doRange = true;
    if length(param.massRange) == 2
        param.massRange = param.massRange(1):param.massRange(end);
    end
end

if ischar(fnames)
    fnames = {fnames};
end

nfiles = size(fnames,2);
ticksz = ceil(nfiles/100);
iupdate = ticksz;
fileSpecs = nan(nfiles,1);
filePars = nan(nfiles,1);
ind = 1;
ppmThreshold = 1; % Not implemented yet

if isempty(dataPath)
    delim = '';
elseif contains(dataPath,'\')
    delim = '\';
else
    delim = '/';
end

for i = 1:nfiles
    if isempty(dataPath)
        ifile = fnames(i);
    else
        ifile = [dataPath delim fnames{i}];
    end

    if i == 1
        % Check file paths
        if regexp(ifile,"_p.h5")
            isProcessed = true;
            avFile = regexprep(ifile,"[/\\]Processed","");
            avFile = regexprep(avFile,"_p","_av");
            noAv = false;
            if ~isfile(avFile)
                avFile = regexprep(avFile,"_av","");
                noAv = true;
            end
            ifFile = regexprep(ifile,"Processed","IF");
            ifFile = regexprep(ifFile,"_p","_IF");
            noIF = false;
            if ~isfile(ifFile)
                noIF = true;
            end
        else
            avFile = ifile;
            ifFile = ifile;
        end
        % ifile
        % avFile
        % ifFile
        % h5disp(ifile)
        % h5disp(avFile)
        % h5disp(ifFile)

        % For the first file, initialize the matrices.
        speci = squeeze(h5read(avFile,'/FullSpectra/TofData'));
        sampleInterval = double(h5readatt(avFile,'/FullSpectra','SampleInterval'));
        startDelay = double(h5readatt(avFile,'/TimingData','StartDelay'));
        tof = startDelay+(0:(length(speci)-1))*sampleInterval*1e9;
        pari = squeeze(h5read(ifFile,'/MassCalib/pars'));
        npar = size(pari,2);

        % Cut away useless data
        if doRange
            if size(pari,1) == 2
                mz1 = ((tof'-pari(2,1))./pari(1,1)).^2;
            elseif size(pari,1) == 3
                mz1 = ((tof'-pari(2,:))./pari(1,:)).^(1./pari(3,:));
            else
                error("Unknown mass calibration")
            end
            iInc = mz1 > (param.massRange(1)-2) & mz1 < (param.massRange(end)+2);
            speci = speci(iInc,:);
            tof = tof(iInc);
        end

        % Initialize matrix for storing spectra
        len = length(speci);
        nspec = size(speci,2);
        ii = ind:(ind+nspec-1);
        specsAll = nan(len,nfiles*nspec);
        specsAll(:,ii) = speci;

        % Mass calibration parameters
        ip = ind:ceil(nspec/npar):(ind+nspec-1);
        pars = nan(size(pari,1),nfiles*nspec);
        pars(:,ip) = pari;

        if doTime
            time = nan(nfiles*nspec,1);
            timeLog = h5read(avFile,'/AcquisitionLog/Log');
            tp = double(timeLog.timestamp/1e7);
            ts = floor(linspace(tp(1),tp(end),nspec+1));
            time(ii) = ts(1:end-1);
            %time(ii) = timeLog.timestamp;
        end
        
        % Additional information that may be used by the algorithm
        %mcalFun = h5readatt(ifile,'/FullSpectra','MassCalibration Function');
        %singleIonSignal = h5readatt(ifile,'/FullSpectra','Single Ion Signal');
        peakShape = h5read(ifFile,'/Peaks/PeakShape');
        widthFit = h5read(ifFile,'/Peaks/PeakWidthFit');
    else
        if isProcessed
            avFile = regexprep(ifile,"[/\]Processed","");
            if ~noAv
                regexprep(avFile,"_p","_av");
            end
            if noIF
                ifFile = ifile;
            else
                ifFile = regexprep(ifile,"Processed","IF");
                ifFile = regexprep(ifFile,"_p","_IF");
            end
        else

        end

        % Read spectra
        speci = squeeze(h5read(ifile,'/FullSpectra/TofData'));
        nspec = size(speci,2);
        ii = ind:(ind+nspec-1);
        % If the number of spectra is bigger than expected
        % expand the data matrices
        if (ind+nspec-1) > size(specsAll,2)
            specsAll = [specsAll, nan(size(specsAll))]; %#ok<AGROW>
            pars = [pars, nan(size(pars))]; %#ok<AGROW>
            if doTime
                time = [time; nan(size(time))]; %#ok<AGROW>
            end
        end
        if doRange
            specsAll(:,ii) = speci(iInc,:);
        else
            specsAll(:,ii) = speci;
        end
        if doTime
            timeLog = h5read(ifile,'/AcquisitionLog/Log');
            tp = double(timeLog.timestamp/1e7);
            ts = floor(linspace(tp(1),tp(end),nspec+1));
            time(ii) = ts(1:end-1);
        end
        % Read mass calibration parameters
        pari = squeeze(h5read(ifile,'/MassCalib/pars'));
        npar = size(pari,2);
        ip = ind:ceil(nspec/npar):(ind+nspec-1);
        pars(:,ip) = pari;
    end
    % Record the number of spectra and parameters in the file.
    fileSpecs(i) = nspec;
    filePars(i) = npar;
    ind = ind+nspec;
    % Update waitbar
    if i == iupdate
        iupdate = iupdate+ticksz;
        if isGUI
            param.guiBar.Value = i/nfiles/2;
        else
            waitbar(i/nfiles,h)
        end
    end
end
% Remove unused cells in datamatrices
irm = all(isnan(specsAll),1);
specsAll(:,irm) = [];
pars(:,irm) = [];
if doTime
    time(irm) = [];
    param.time = time;
end
parsFinal = fillmissing(pars,'previous',2);
if isGUI
    param.guiBar.Value = 1/2;
    param.guiBar.Message = "Preparing mass axes...";
else
    waitbar(1,h,"Preparing mass axes...")
end

% Create the mass axes from calibration parameters
if size(pars,1) == 2
    mzAll = ((tof'-parsFinal(2,:))./parsFinal(1,:)).^2;
elseif size(pars,1) == 3
    mzAll = ((tof'-parsFinal(2,:))./parsFinal(1,:)).^(1./parsFinal(3,:));
else
    error("Unknown mass calibration")
end

if isGUI
    param.guiBar.Value = 2/3;
else
    close(h);
end

% Read parameters from imported data
param.peakShape.dat = peakShape';
param.peakShape.dat(1:3,2) = 0;
param.peakShape.dat(end-2:end,2) = 0;% Force slope to zero
Rx = (1:(length(widthFit)-1))';
Rpar = polyfit(Rx,Rx./widthFit(2:end),1);
param.R = @(m) m./(Rpar(1)*m+Rpar(2));
param.Rpar = Rpar;

% Only change this if you really know what you're doing. The value 0
% results in the algorithm automatically determining the baseline.
param.baseline = 0;
% Only change this if you really know what you're doing.
param.weights = 1;

end

