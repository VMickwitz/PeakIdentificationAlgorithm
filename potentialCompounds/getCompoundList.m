function [mass,name,nel] = getCompoundList(elements,ll,ul,charge)
%GETCOMPOUNDLIST Summary of this function goes here
%   Detailed explanation goes here
% 

n = length(elements);
limits = [ll(:) ul(:)];
elmass = getMass(elements);
[elmass,isort] = sort(elmass);
limits = limits(isort,:);

nel = getAllCombinations(limits);
nel(:,isort) = nel;
mass = nel*elmass;
[mass,isort] = sort(mass);
nel = nel(isort,:);

if nargin < 4
    charge = "";
end

electronMass = 0.000548579909;
if strcmp(charge,"-")
    mass = mass+electronMass;
elseif strcmp(charge,"+")
    mass = mass-electronMass;
end

name = strings(length(mass),1);
for i = 1:length(mass)
    for j = 1:n
        if nel(i,j) > 0
            name(i) = sprintf("%s%s%i",name(i),elements(j),nel(i,j));
        end
    end
    name(i) = sprintf("%s%s",name(i),charge);
end
    

end


function [masses] = getMass(names)
masses = nan(length(names),1);
for i = 1:length(names)
    m_i = getIsoMass(names(i),0.1);
    masses(i) = m_i(1);
end

end