function showAssignment(fit,M)
%SHOWASSIGNMENT Shows process from preliminary to final fit.

compPar.verbate = 1;
param = fit.param;
if isfield(param,'compF')
    if isstruct(param.compF)
        % if compF is a list, reformat into function:
        param.compF = @(x,par) findComp(x,param.compF,par);
    % elseif isa(param.compF,'function_handle')
    %     param.compF = @(x,par) param.compF(x,par);
    end
else
    error("No compF provided.")
end

i = find(param.massRange == M,1,'first');
i_M = round(fit.mz)==M;

fit_i = fit;
fit_i.param.compF = param.compF;
fit_i.param.W = param.W(i);
fit_i.param.A = param.A(i);
fit_i.param.baseline = fit.param.baseline(i_M,:)./param.W(i);
fit_i.mz = fit.mz(i_M);
fit_i.specs = (fit.specs(i_M,:)-fit.isoSpecs(i_M,:))./param.W(i);
isoY = fit.isoSpecs(i_M,:);
fit_i.isoSpecs = isoY./param.W(i);
fit_i.param = rmfield(fit_i.param,'BL');


fitAndAssign(fit_i,1,M,compPar);

end

