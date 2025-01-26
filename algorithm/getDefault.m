function [value] = getDefault(fieldName,fit)
%GETDEFAULT Returns default value of parameter specified in fieldName.

switch fieldName
    case "A"
        value = 1;
    case "baseline"
        if nargin > 1
            value = zeros(size(fit.specs));
        else
            error("This field requires second input")
        end
    case "compF"
        value = @(x) find_composition(x);
    case "dmcal"
        value = 0;
    case "doRm"
        value = false;
    case "massRange"
        if nargin > 1
            value = sort(unique(round(mz)));
        else
            error("This field requires second input")
        end
    case "peakRange"
        value = 0:10;
    case "resFactor"
        value = 1.2;
    case "Rtol"
        value = 0;
    case "tol"
        value = 0;
    case "weights"
        if nargin > 1
            value = ones(1,size(specs,2));
        else
            error("This field requires second input")
        end
    otherwise
        value = -1;
        % Used for error processing in checkParam.
end


