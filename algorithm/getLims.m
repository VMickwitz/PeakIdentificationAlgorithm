function [Alims,numbers,chis,irm] = getLims(numbers,chis)
%Calculates limits of best A values based on numbers of peaks and chi.
%   Final result has to be multiplied by factor ln(length(x)*size(y,2)) for
%   consistency with other algorithm functions.

numbers_in = numbers;
Alims = nan(length(numbers),2);
Alims(1,1) = 0;
i = 1;
while i < length(numbers)
    if chis(i)-chis(i+1) > 0
        % Change the upper limit of current index, and lower limit of
        % coming index.
        Alims(i,2) = (numbers(i+1)-numbers(i))./(chis(i)-chis(i+1));
        Alims(i+1,1) = Alims(i,2);
        if i > 0 && Alims(i,2) <= Alims(i,1)
            % If new upper limit is lower than lower limit, the current
            % idex is never chosen. Remove it.
            numbers(i) = [];
            Alims(i,:) = [];
            Alims(i,1) = Alims(i-1,2);
            chis(i) = [];
            i = i-1;    % Move back to previous index.
        else
            i = i+1;    % Move one index forward.
        end
    else
        % If following index is never an improvement over the current one.
        % Remove it.
        numbers(i+1) = [];
        Alims(i+1,:) = [];
        chis(i+1) = [];
        % Stay on current index.
    end
end

% Make final line of lims:
Alims(end,2) = inf;

[~,irm] = setdiff(numbers_in,numbers);
end

