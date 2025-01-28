
%% Clear other variables. Can be skipped.
% Lines marked with "%Edit" should be changed by the user.
% Lines marked with "%Change" can be changed by an experienced user.

% If you improve any functions related to the algorithm, please consider
% uploading it to github or sending it to valter.mickwitz@helsinki.fi for
% the benefit of other users.

clear;

%% Specify folder for data
% Currently all files in the folder are used for the algorithm
% It is recommended to copy the files you want to use as reference for
% peaklist generation to a new folder, and specify the path here. This will
% ensure your original data files are not lost due to bugs.

% This will open a gui asking you to provide the files with your data.
% Make sure your files contain mass calibration and peak shape data.
[files, dataPath] = uigetfile(["*.h5" "*.hdf"],'MultiSelect','on');

% .mat file where the results will be saved. If a file already exists it
% will be overwritten. Be careful.
saveName = "autofit_test.mat"; % Edit

addAppPaths(); % Adds relevant paths to matlab path.

%% Read information from files
% If reading in a lot of data matlab may run out of memory here.
% This may get improved in the future. For now, try to restrict the data
% used for peaklist generation if this part gives an error.

[specsAll,mzAll,param] = importDataset(dataPath,files);

%% Define parameters for algorithm, these can all be tweaked.

% massRange determines whitch integer masses will be analyzed
param.massRange = 300:310; % Edit

% lims determine the window within which peaks are searched around an
% integer mass, eg. [-0.2 0.3] means the algorithm looks for peaks at
% integer mass 355 between 354.8 and 355.3.
param.lims = [-0.3 0.4]; % Edit

% peakRange determines the minimum and maximum number of peaks that may be
% present at an integer mass. Higher max values result in slower runs, but
% max should be higher than what you would expect to find at any mass.
% Recommended to keep lower limit at 0.
param.peakRange = [0 11]; % Edit

% compF is the list of potential compounds to use.
param.compF = importdata("potentialListCIMS.mat"); % Edit

% Limits for the what the average number of peaks per unit mass may be.
% This restricts the values of parameter A during the algorithm running.
param.avgPeakLims = [3 10]; % Change

%% Create the fit structure. Interpolate and average the imported data.

fit = struct;
[fit.mz,specs] = interpmz(specsAll,mzAll,param);
fit.specs = sum(specs,2);
fit.param = param;

%% Plot spectrum to check the imported data looks ok.

figure;
plot(fit.mz,sum(specs,2),'k-','LineWidth',2)
xlim([param.massRange(1), param.massRange(end)])

%% Uncomment to check peakShape

% figure;
% plot(param.peakShape.dat(1,:),param.peakShape.dat(2,:),'k-','LineWidth',2)
%% Run the algorithm
verb = true;    % Change (boolean, determines amount of output)
doPre = true;   % Change (boolean, determines if the algorithm should run 
                % the preliminary fitting or not).
fit = autoFit(fit,saveName,verb,doPre);
save(saveName,"fit")
%% Plots results for individual masses and displays the potential compounds

% This section can be used to check results at individual masses before 
% exporting the peaklist.
M = 340; % Change
list = param.compF;

if M < fit.param.massRange(1) || M > fit.param.massRange(end)
    warning("Mass outside of analyzed range. No plot generated.")
else
    plotMass(fit,M,'sparseLegend',true);
    ind = round(fit.PeakList) == M;
    iList = round(list.mass) == M;
    locatedCompounds = [string(fit.PeakList(ind)) fit.CompNams(ind)] %#ok<NOPTS>
    potentialCompounds = [string(list.mass(iList)) list.names(iList)] %#ok<NOPTS>
end

%% Export generated peaklist

% Set filepath manually or leave empty to be prompted by gui
filepath = strrep(saveName,".mat","")+"_list.txt"; %Change
%filepath = [];

% Currently the only available export format is "tofware".
% If you want a different format you can edit the exportPeakList function
% or contact valter.mickwitz@helsinki.fi
format = "tofware";
exportPeakList(fit,filepath,format);
% The generated .txt file can be directly imported to the software
% specified by the format input for further analysis.


