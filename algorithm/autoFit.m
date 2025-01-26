function [fitStruct] = autoFit(fit,saveNam,verb,doPre)

if nargin <2 || isempty(saveNam)
    doSave = false;
else
    doSave = true;
end
if nargin < 3 || isempty(verb)
    verb = false;
end
if nargin < 4 || isempty(doPre)
    % Used mainly for testing. Can load a previous perliminary fit insted
    % of refitting.
    doPre = true;
end

% Get preliminary fit:
if doPre
    param = fit.param;
    fit = rmfield(fit,'param');
    fitPre = fitForA(fit,param);
    if doSave
        save(saveNam,"fitPre")
    end
else
    if doSave
        fit = importdata(saveNam);
    end
    fitPre = fit;
end

if verb
    fitPre %#ok<NOPRT>
    fitPre.param
    fitPre.fitPrelim
end

% Determine A
fitPre = getA(fitPre);

fitStruct = constructPeakList(fitPre);
if verb
    fprintf("Final fit structure:")
    fitStruct %#ok<NOPRT>
end
if doSave
    save(saveNam,"fitStruct")
end

end

