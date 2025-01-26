function [N] = getAllCombinations(limits)
%GETALLCOMBIANTIONS returns all combinations of integer element values
%between (including) limits.
% First column of limits is lower limit, second is upper. Integer limits
% assumed

limits = int32(limits);
if any(limits<0)
    error("No negative limits allowed.")
elseif any(limits(:,1)>limits(:,2))
    limits
    error("Check limits.")
end

n = length(limits);
delta = limits(:,2)-limits(:,1)+1;

ntot = prod(delta);
N = zeros(ntot,n);
Ni = ntot/delta(1);
N(:,1) = repmat((limits(1,1):limits(1,2))',[Ni,1]);
Nj = delta(1);
for i = 2:n
    Ni = Ni/delta(i);
    vec = repmat(limits(i,1):limits(i,2),[Nj,Ni]);
    N(:,i) = vec(:);
    Nj = Nj*delta(i);
end



end

