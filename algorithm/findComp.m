function [names,masses,ind] = findComp(mass,list,compPar)
% Find composition close to input mass from sorted list.
%   Inputs:
%       mass - scalar value for where to look for compounds.
%       list - structure with fileds:
%           mass: masses of compositions.
%           names: names of compositions.
%       compPar - structure with parameters for composition
%           tol: maximum difference in mass from input mass.

% Tolerance for mass
tol = 0.01;
if nargin > 2
    if isfield(compPar,'tol')
        tol = compPar.tol;
    end
end

% Find indices for masses within tol:
lims = mass + [-tol tol];
ind1 = binarySearch(list.mass,lims(1));
ind2 = binarySearch(list.mass,lims(2),ind1)-1;

% Select correct indices to return:
masses = list.mass(ind1:ind2);
names = list.names(ind1:ind2);

if nargout > 2
    % Optionally return indices
    ind = [ind1, ind2];
end

end

