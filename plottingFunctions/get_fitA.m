function new_fit = get_fitA(fit,A,fieldName)

% If a limit is equal to A, the lower number of peaks is selected.

doField = false;
if nargin > 2
    doField = true;
end

lims = fit.Alims;

n = size(lims,3);

peaks_all = fit.peaks;
h_all = fit.H;
numbers_all = fit.numbers;
nmin = fit.numbers(1,1);
doBLs = false;

if isfield(fit,'BLs')
    % This has not been tested
    BLs_all = fit.BLs;
    doBLs = true;
end

if isfield(fit,'mapA') && fit.mapA
    lims = mappingFunction(permute(lims,[3,1,2]),fit.param.W);
    lims = permute(lims,[2 3 1]);
    %b = 2.7;
    %lims = permute(permute(lims,[3,1,2])./exp(log(fit.param.W) - b*(1 + exp(-fit.param.W.^0.5))),[2,3,1]);
end

peaks = nan(n,size(peaks_all,2));
h_fit = nan(size(h_all,1),size(h_all,2),n);
numbers = nan(n,1);
BLs = nan(n,1);
if doField
    value = nan(n,1);
    value_all = fit.(fieldName);
end

if length(A) == 1
    A = A*ones(n,1);
end


for i = 1:n
    if A == 0
        ind = 1;
    else
        ind = squeeze(lims(:,1,i) < A(i) & lims(:,2,i) >= A(i));
    end
    peaks(i,:) = peaks_all(ind,:,i);
    h_fit(:,:,i) = h_all(:,:,ind,i);
    numbers(i) = numbers_all(ind,i);
    if doField
        % Indexing not consistent. Careful.
        value(i,:) = value_all(ind,i);
        %value(i) = value_all(numbers(i)+1-nmin,i);
    end
    if doBLs
        BLs(i) = BLs_all(ind,i);
    end
end

new_fit = fit;
new_fit.lims = lims;
new_fit.peaks = peaks;
new_fit.H = h_fit;
new_fit.numbers = numbers;
new_fit.W = fit.param.W;
new_fit.metaFit = fit;
new_fit.BLs = BLs;
if doField
    new_fit.(fieldName) = value;
end
end

