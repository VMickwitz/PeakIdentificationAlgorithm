function [mass,elnum,name] = buildComp(elnam,seeds,clusters,nmax,mlim,charge)
%BUILDCOMP Summary of this function goes here
%   Detailed explanation goes here

if nargin < 5
    charge = "";
elseif ~isstring(charge)
    charge = string(charge);
end

load("massStruct.mat","massStruct")
nel = length(elnam);
elmass = nan(nel,1);
for i = 1:nel
    if elnam(i) == "X"
        % Ficticious element used for algorithm testing.
        elmass(i) = 1;
    else
        elmass(i) = massStruct.(elnam(i));
    end
end

% Calculate cluster and seed masses
cmass = clusters*elmass;
smass = seeds*elmass;
if any(smass<=0) || any(cmass<=0)
    error("No negative cluster masses allowed.")
end
cmin = min(cmass);
cmax = max(cmass);
smin = min(smass);
smax = max(smass);
% Make some upper and lower limits for total cluster number.
maxn = max(floor((mlim(2)-smin)/(cmin)),0);
minn = max(ceil((mlim(1)-smax)/cmax),0);

% Make all possible combinations of cluster numbers.
if length(nmax)==1
    % If only one cluster allowed, just add 0 to max.
    nclusts = (0:nmax)';
else
    if ~(length(nmax)==size(clusters,1))
        error("Maximum cluster number inconsistent with number of clusters.")
    end
    nmax = nmax(:);
    limits = [zeros(length(nmax),1), nmax];
    % Get all combinations of numbers within limits:
    nclusts = getAllCombinations(limits);
end

% Remove too big/small masses
ntot = sum(nclusts,2);
nclusts(ntot>maxn|ntot<minn,:) = [];

% Consturct the element number matrix:
elnums = nclusts*clusters;
seedmat = repmat(seeds,1,size(elnums,1))';
elnums = repmat(elnums,size(seeds,1),1);
elnum = elnums+reshape(seedmat,nel,size(elnums,1))';
elnum = unique(elnum,"rows");

% Calculate mass of compounds:
mass = elnum*elmass;

% Account for charge
electronMass = 0.000548579909; %Da
if strcmp(charge,"-")
    mass = mass+electronMass;
elseif strcmp(charge,"+")
    mass = mass-electronMass;
elseif ~strcmp(charge,"")
    error("Failed to read charge.")
end

% Remove masses with negative element counts or outside mass range.
irm = mass<mlim(1) | mass>mlim(2) | any(elnum<0,2);
mass(irm) = [];
elnum(irm,:) = [];
ncomp = size(elnum,1);

% Make the names:
pattern = strtrim(string(repmat('%s%i ',1,nel)));
ivec = 1:2:2*nel;
strcell = cell(ncomp,2*nel);
namcell = repmat(cellstr(elnam),ncomp,1);
[strcell{:,ivec}] = deal(namcell{:});
[strcell{:,2*nel+1}] = deal(charge);
numcell = num2cell(elnum(:));
[strcell{:,ivec+1}] = deal(numcell{:});
strcell = strcell';
name = strsplit(sprintf(pattern+"%s"+"/",strcell{:}),"/");
name = name(1:end-1)';

% Remove zeroes from names.
name = regexprep(name,"[A-Z][a-z]?[0][ ]?","");
name = regexprep(name," ["+charge+"]",charge);
% Sort by mass:
[mass,isort] = sort(mass);
elnum = elnum(isort,:);
name = name(isort);

end