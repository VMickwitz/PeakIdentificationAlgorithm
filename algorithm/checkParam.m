function [value] = checkParam(fit,fieldName)
%CHECKPARAM If value not in param, returns default,

if isfield(fit.param,fieldName)
    value = fit.param.(fieldName);
else
    value = getDefault(fieldName,fit);
end
end

