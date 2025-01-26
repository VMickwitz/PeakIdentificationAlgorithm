function [sorted1,sorted2,ind1,ind2] = custom_sort(arr1,arr2,maxDiff)
% Attempts to pair the values in two different sized arrays. Adds nan
% values to fill in gaps. Not made for large arrays, very basic, has some bugs.
%   arr1 = sorted1(ind1)
%   Detailed explanation goes here

arr1(isnan(arr1)) = [];
arr2(isnan(arr2)) = [];

if nargin < 3
    maxDiff = 0.01;
end

if length(arr1) > length(arr2)
    flip = false;
    long = arr1;
    short = arr2;
else
    flip = true;
    long = arr2;
    short = arr1;
end

[long,iL] = sort(long,'ascend');
[short,iS] = sort(short,'ascend');
[~,iL] = sort(iL);
[~,iS] = sort(iS);
newlong = nan(length(long)+length(short),1);
newshort = newlong;
indL = nan(length(long),1);
indS = nan(length(short),1);

k = 1;
indOld = 0;

for i = 1:length(short)
    % Check that the two compared values are each other nearest candidate
    [~,ind1] = min((long-short(i)).^2);
    [~,ind2] = min((short-long(ind1)).^2);
    % Check if some values in long were skipped, and add them
    if ind1 > indOld+1
        newshort(k:k+ind1-indOld-2) = nan;
        newlong(k:k+ind1-indOld-2) = long(indOld+1:ind1-1);
        indL(indOld+1:ind1-1) = k:k+ind1-indOld-2;
        k = k+ind1-indOld-1;
        indOld = ind1-1;
%         fprintf("1\n")
    end
    if ind2 == i
%         % Check if some values in long were skipped, and add them
%         if ind1 > indOld+1
%             newshort(k:k+ind1-indOld-2) = nan;
%             newlong(k:k+ind1-indOld-2) = long(indOld+1:ind1-1);
%             indL(indOld+1:ind1-1) = k:k+ind1-indOld-2;
%             k = k+ind1-indOld-1;
%         end
        % Check that difference is below maxDiff
        if abs(short(i)-long(ind1)) < maxDiff
            % Pair values
            newshort(k) = short(i);
            newlong(k) = long(ind1);
            indS(i) = k;
            indL(ind1) = k;
            k = k+1;
            indOld = ind1;
%             fprintf("2\n")
        % Add entries in correct order if not paired
        elseif short(i)-long(ind1) < 0
            newshort(k) = short(i);
            newlong(k) = nan;
            indS(i) = k;
            newshort(k+1) = nan;
            newlong(k+1) = long(ind1);
            indL(ind1) = k+1;
            indOld = ind1;
            k = k+2;
%             fprintf("3\n")
        else
            newshort(k) = nan;
            newlong(k) = long(ind1);
            indL(ind1) = k;
            newshort(k+1) = short(i);
            newlong(k+1) = nan;
            indS(i) = k+1;
            indOld = ind1;
            k = k+2;
%             fprintf("4\n")
        end
    else
        newshort(k) = short(i);
        newlong(k) = nan;
        indS(i) = k;
        k = k+1;
%         fprintf("5\n")
    end
%     newshort
%     newlong
end

if indOld < length(long)
    newlong(k:k+length(long)-indOld-1) = long(indOld+1:end);
    indL(indOld+1:end) = k:k+length(long)-indOld-1;
    k = k+length(long)-indOld;
%     fprintf("6\n")
end

%remove nan values at the end.
newlong(k:end) = [];
newshort(k:end) = [];

% check array
check = true;
while check
    check = false;
    for i = 1:length(newlong)-1
        
    end
    
if ~flip
    sorted1 = newlong;
    sorted2 = newshort;
    ind1 = indL(iL);
    ind2 = indS(iS);
else
    sorted1 = newshort;
    sorted2 = newlong;
    ind1 = indS(iS);
    ind2 = indL(iL);
end
    
end

