function [iso] = getIsoMass(comp,pmin)
%GETISOMASS Returns isotope masses and relative abundance for given
%compound.

% It's faster to save the information for every number of atoms of each
% element, and look it up when needed, than calculating.

% Inputs:
%   comp - composition for which to get the isotope masses and abundances
%   pmin - (optional) minimum relative abundance of output isotopes

% pmin = 1 just returns the composition mass consisting of the most common
% isotopes.

rmP = true;
if nargin < 2
    pmin = 1e-7;
end

load("isoStruct.mat","isoStruct")


if length(comp) > 1
    error("Only one composition can be handled at a time.")
elseif isa(comp,"cell")
    comp = string(comp{1});
end

charge = regexp(comp,"([+-]$)",'tokens'); % Charge must be at the end of the input string.
matches = regexp(comp,"([A-Z][a-z]?)(\d*)",'tokens');

electronMass = 0.000548579909; %Da

for i = 1:length(matches)
    elinfo = matches{i};
    if strcmp(elinfo(2),"")
        nel = 1;
    else
        nel = str2double(elinfo(2));
    end

    if strcmp(elinfo(1),"X")
        % For testing purposes. X has one isotope only, with mass 1u.
        m = nel;
        p = 1;
    elseif nel == 1
        isoDat = isoStruct.(elinfo(1)); 
        p = isoDat(:,3);
        m = isoDat(:,2);
    elseif isfield(isoStruct,elinfo(1)+elinfo(2))
        iso_i = isoStruct.(elinfo(1)+elinfo(2));
        m = iso_i(:,1);
        p = iso_i(:,2);
    else
        p = 1;
        isoDat = isoStruct.(elinfo(1));
        p0 = isoDat(:,3)'; 
        if isscalar(p0)
            m = nel*isoDat;
        else
            niso = zeros(size(p0));
            for j = 1:nel
                % Iteratively add elements, and update probability of each
                % chain.
                % p and n change size every iteration. Fixable?
                [p,niso] = addEl(p,niso,p0);
            end
            for k = 1:length(p)
                % Account for multinomial coefficient
                npart = nchooseks(niso(k,:));
                p(k) = p(k).*npart;
                m = sum(niso.*isoDat(:,2)',2);
            end
        end
    end
    if i == 1
        % On first iteration create the array.
        if pmin == 1
            iso = [m(1), p(1)];
        else
            iso = [m, p];
        end
    elseif pmin == 1
        [iso] = mergeElements(iso,p(1),m(1));
    else
        % Merge this element count to other previous elements.
        [iso] = mergeElements(iso,p,m);
        if rmP
            iso(iso(:,2) < pmin,:) = [];
        end
    end
end

[iso(:,2),is] = sort(iso(:,2),'descend');
iso(:,1) = iso(is,1);

if isempty(charge)
    %warning("No charge")
elseif strcmp(charge{1},"-")
    iso(:,1) = iso(:,1) + electronMass;
elseif strcmp(charge{1},"+")
    iso(:,1) = iso(:,1) - electronMass;
else
    warning("Charge detemination failed.")
end

end 

function [p,n] = addEl(p,n,p0)
% p is column vector of previous probabilities
% p0 is row vector of isotope relative abundancies
np = length(p0);
nlen = length(p);
p = p(:).*p0;
%n = repmat(n,[np,1]);
n = repmat(reshape(n',[],1),[1,np]);
n = n + repmat(diag(ones(np,1)),[nlen,1]);
[n,in] = unique(reshape(n,[np,np*nlen])','rows');
p = p(in);
end

function [npartitions] = nchooseks(ns)
% The number of ways sum(ns) elements can be divided into length(ns)
% groups, with group(i) containing ns(i) elements.
nsum = cumsum(ns);
npartitions = 1;
for ind = 2:length(ns)
    if ns(ind)
        npartitions = npartitions*nchoosek(nsum(ind),ns(ind));
    end
end

end

function [iso] = mergeElements(iso,p,m)
m = iso(:,1)+m';
p = iso(:,2).*p';
iso = [m(:), p(:)];
end