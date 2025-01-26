function [idx] = binarySearch(A,value,i1,i2)
%Finds first index idx stisfying A(idx) >= value, for sorted array A.

% Initialize i1 and i2 if required
if nargin < 3
    i1 = 1;
    i2 = length(A);
elseif nargin == 3
    i2 = length(A);
end

% Iteratively improve i1 (lower limit)
% and i2 (upper limit)
while i1 < i2
    idx = ceil((i1+i2+1)/2)-1;
    if A(idx) < value
        i1 = idx+1;
    else
        i2 = idx;
    end
end

% Check that i1 and i2 converged to same value.
if i1 == i2
    idx = i1;
end

% If value is out of range, return -1.
if A(idx) < value
    idx = -1;
end

end

