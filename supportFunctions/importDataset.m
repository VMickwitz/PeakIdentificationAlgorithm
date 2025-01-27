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

if isfield(param,'guiBar')
    isGUI = true;
    param.guiBar.Message = "Importing Spectra...";
    param.guiBar.Value = 0;
else
    h = waitbar(0,"Importing spectra...");
end
nfiles = length(fnames);
ticksz = ceil(nfiles/100);
iupdate = ticksz;
fileSpecs = nan(nfiles,1);
filePars = nan(nfiles,1);
ind = 1;
ppmThreshold = 1; % Not implemented yet

if contains(dataPath,'\')
    delim = '\';
else
    delim = '/';
end

for i = 1:nfiles
    ifile = [dataPath delim fnames{i}];
    if i == 1
        % For the first file, initialize the matrices.
        speci = squeeze(h5read(ifile,'/FullSpectra/TofData'));
        sampleInterval = double(h5readatt(ifile,'/FullSpectra','SampleInterval'));
        startDelay = double(h5readatt(ifile,'/TimingData','StartDelay'));

        % Initialize matrix for storing spectra
        len = length(speci);
        nspec = size(speci,2);
        ii = ind:(ind+nspec-1);
        specsAll = nan(len,nfiles*nspec);
        specsAll(:,ii) = speci;

        % Mass calibration parameters
        pari = squeeze(h5read(ifile,'/MassCalib/pars'));
        npar = size(pari,2);
        ip = ind:ceil(nspec/npar):(ind+nspec-1);
        pars = nan(size(pari,1),nfiles*nspec);
        pars(:,ip) = pari;
        
        % Additional information that may be used by the algorithm
        %mcalFun = h5readatt(ifile,'/FullSpectra','MassCalibration Function');
        %singleIonSignal = h5readatt(ifile,'/FullSpectra','Single Ion Signal');
        peakShape = h5read(ifile,'/Peaks/PeakShape');
        widthFit = h5read(ifile,'/Peaks/PeakWidthFit');
    else
        
        % Read spectra
        speci = squeeze(h5read(ifile,'/FullSpectra/TofData'));
        nspec = size(speci,2);
        ii = ind:(ind+nspec-1);
        % If the number of spectra is bigger than expected
        % expand the data matrices
        if (ind+nspec-1) > size(specsAll,2)
            specsAll = [specsAll, nan(size(specsAll))]; %#ok<AGROW>
            pars = [pars, nan(size(pars))]; %#ok<AGROW>
        end
        specsAll(:,ii) = speci;
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
parsFinal = fillmissing(pars,'previous',2);
if isGUI
    param.guiBar.Value = 1/2;
    param.guiBar.Message = "Preparing mass axes...";
else
    waitbar(1,h,"Preparing mass axes...")
end

% Create the mass axes from calibration parameters
tof = startDelay+(0:len-1)*sampleInterval*1e9;
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
if isfield(param,'massRange') && length(param.massRange) == 2
    param.massRange = param.massRange(1):param.massRange(end);
end

% Only change this if you really know what you're doing. The value 0
% results in the algorithm automatically determining the baseline.
param.baseline = 0;
% Only change this if you really know what you're doing.
param.weights = 1;

end

