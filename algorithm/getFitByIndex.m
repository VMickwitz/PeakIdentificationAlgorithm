function [n, peaks, H, BL] = getFitByIndex(fitStruct,ind,M)
% Returns the indexed fit in order of goodness.

mz = fitStruct.mz(round(fitStruct.mz)==M);
l = length(mz);
s = size(fitStruct.specs,2);
i = find(fitStruct.param.massRange==M,1,'first');

% fprintf("by Index\n")
% fitStruct.param.A
% fitStruct.chi(:,i)
% fitStruct.numbers(:,i)
% fitStruct.param.peakRange

%scores = fitStruct.param.A*fitStruct.chi(:,i)+fitStruct.numbers(:,i).*log(s*l);

scores = fitStruct.param.A*fitStruct.fitPrelim.chi(:,i)+fitStruct.fitPrelim.numbers(:,i).*log(s*l);

% if isempty(scores)
%     M
%     i
%     ind
%     fitStruct.chi
%     fitStruct.param.peakRange'
%     s
%     l
%     fitStruct.param.A
% end

%[~,iSort] = sort(scores,1,'ascend');
[~,ind] = min(scores);
n = fitStruct.fitPrelim.numbers(ind,i);
ind = find(fitStruct.fitPrelim.numbers(:,i) == n,1,'first');

peaks = fitStruct.fitPrelim.peaks(ind,:,i);
H = fitStruct.fitPrelim.H(:,:,ind,i);

peaks(isnan(peaks)) = [];
H(isnan(H)) = [];

if isfield(fitStruct.fitPrelim,'BLs')
    BL = fitStruct.fitPrelim.BLs(ind,i);
end

end

