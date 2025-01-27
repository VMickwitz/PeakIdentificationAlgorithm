function [mzout,specsout] = interpmz(specs,mz,param)
%INTERPMZ Interpolates data to common massaxis. Also trims all values
%outside range specified in parameters.

lims = param.lims;
range = param.massRange;
mzout = mz(:,1);
doInterp = true;
if size(mz,2) == 1
    doInterp = false;
end

isGUI = false;
if isfield(param,'guiBar')
    isGUI = true;
    param.guiBar.Message = "Interpolating spectra to common axis...";
else
    h = waitbar(0,"Interpolating spectra to common axis...");
end

% Remove all useless data:
include = false(size(mzout));
for m = range
    %mzout<(m+lims(2)) & mzout>(m+lims(1))
    include(mzout<(m+lims(2)) & mzout>(m+lims(1))) = true;
end

nspecs = size(specs,2);
stepsz = floor(nspecs/100);
iupdate = stepsz;

mzout = mzout(include);

if ~doInterp
    specsout = specs(include,:);
    return;
end

specsout = nan(sum(include),size(specs,2));
specsout(:,1) = specs(include,1);

for i = 2:nspecs
    i1 = find(mz(:,i) > (range(1)-1),1,'first');
    i2 = find(mz(:,i) < (range(end)+1),1,'last');
    mzuse = mz(i1:i2,i);
    specuse = specs(i1:i2,i);
    specsout(:,i) = spline(mzuse,specuse,mzout);
    if i == iupdate
        iupdate = iupdate+stepsz;
        if isGUI
            param.guiBar.Value = 2/3 + i/nspecs/3;
        else
            waitbar(i/nspecs)
        end
    end
end

if ~isGUI
    close(h)
end
end
