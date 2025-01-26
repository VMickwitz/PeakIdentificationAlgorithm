function [masses, compNam, doAssign, AltMasses, AltComp] = assignPeaks(mz,specs,compList,param,compPar)
%ASSIGN_PEAKS Assigns a formula to the first peak in compList
%   This is an old function so it is probably needlessly complicated and
%   contains some rudimentary lines from earlier versions of the algorithm.
%   Inputs:
%       mz - mass axis at unit mass
%       specs - spectra at unit mass
%       compList - list of masses of free peaks
%       param - parameter structure
%       compPar - parameters specifically for formula assignment

% Initilaization
masses = [];
compNam = [];
AltMasses = [];
AltComp = {};
param1 = param;
useCompF = false;
assignDiff = 0.005;
verbate = false;

% Read relevant parameters
if isfield(param,'peakList')
    % List of all peaks
    param1.peakList = [compList param1.peakList];
else
    param1.peakList = compList;
end
if isfield(param,'compF')
    % Function used to get potential compounds
    useCompF = true;
    assert(isa(param.compF,'function_handle'),"Given compF not a function handle.")
else
    error("No compF provided.")
end
if isfield(compPar,'assignDiff')
    % Maximum allowable difference in mass between free peak and assigned
    % formula
    assignDiff = compPar.assignDiff;
end
if isfield(compPar,'verbate')
    % Write outputs
    verbate = compPar.verbate;
end

doAssign = true;

% Read the mass to be assigned a formula.
imax = 1; % Index of free peaks that gets a formula
comps3 = compList; 
if isempty(comps3)
    % If no free peaks: skip
    doAssign = false;
    return;
else
    comps3(imax) = [];% Free peaks that will not get assigned a formula
end

% Get names and masses of potential formulas
if useCompF
    [NamS,compMass] = param.compF(compList(imax),compPar);
    compDif = compMass-compList(imax);
end

if verbate
    compList(imax)
    NamS
    compDif
end


if isempty(compMass)
    % If no potential formulas: skip
    doAssign = false;
else
    % Compare difference in mass to the maximum allowed difference
    NamSimp = NamS';
    Dif = abs(compDif);
    if verbate
        assignDiff
    end
    Inc = Dif < assignDiff;
    % Add formulas outside interval to list of alternative compounds
    AltComp = NamSimp(~Inc);
    AltMasses = compMass(~Inc);
    compNams = NamSimp(Inc);
    compMass = compMass(Inc);


    if isempty(compMass)
        % If no formulas within interval: skip
        doAssign = false;
%    elseif ~isempty(param.peakList) && any(param.peakList==compMass(1))
%         doAssign = false;
    elseif (sum(Inc)==1) && (~any(masses==compMass(1)))
        % If only one formula in interval: assign it.
        masses = [masses compMass(1)];
        compNam = [compNam; string(compNams)];
        if isempty(AltMasses)
            AltMasses = nan;
            AltComp = "none";
        end
    elseif (sum(Inc)>1)% && (~any(masses==compMass(1)))
        % If several potential formulas: compare their effects on fit
        evalF = nan(length(compMass),1);
        for i = 1:length(compMass)
            mass = compMass(i);
            param1.peakList = [param.peakList mass];
            param1.method = 3;
            param1.tol = assignDiff; % Allows free peaks to move within these limits
            if verbate
                [~,fvals,~,~,par] = fitPeaks(mz,specs,comps3,param1);
                mass
                fprintf("%s: fval = %5.4e\n",compNams(i),fvals(1))
            else
                [~,fvals,~,~,~] = fitPeaks(mz,specs,comps3,param1);
            end
            evalF(i) = fvals(1);
            
        end
        [~,iMin] = min(evalF);

        if ~any(masses==compMass(iMin))
            masses = [masses compMass(iMin)];
            compNam = [compNam; compNams(iMin)];
            AltMasses = [compMass(2:end); AltMasses];
            %AltComp = string(compNams{1,2:end});
            compNams2 = compNams;
            compNams2(iMin) = [];
            AltComp = [compNams2, AltComp];
            %AltComp = string(ones(size(compNams2,2),1))
            %for i = 1:size(compNams2,2)
            %    AltComp(i,1) = string(compNams2{1,i})
            %end
            %AltComp = string(compNams2{1,1:end});
            %AltComp = string(AltComp);
        end
    else
        if verbate
            fprintf("No mass in range")
        end
        AltComp = [compNams, AltComp];
        AltMasses = [compMass; AltMasses];
    end

    if ~isempty(masses) && any(param.peakList == masses)
        % Remove compound if too close to already assigned compound, and try again.
        [masses, compNam, doAssign, AltMasses, AltComp] = assignPeaks(mz,specs,comps3,param,compPar);
    end
end
    
end

