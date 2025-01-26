function [isoYs,isoH,Pisos,Hisos,parent] = getIsoYs(fit,compList,M,approx)
%GETISOYS Returns estimated spectra of less common isotopes
compNams = compList.compNams;
peakList = compList.peakList;
useH = false;
if isfield(compList,'Hgen')
    useH = true;
    Hgen = compList.Hgen;
    approx = false;
elseif ~isfield(fit,'specs')
    error("Input fit structure requires specs field when no Hgen given.")
end

if nargin < 4
    % If approx is true function attempts to approximate impact of isotopes
    % on previous peaks.
    approx = false;
end

i0 = round(fit.mz) == M;
x0 = fit.mz(i0);
isoYs = zeros(length(x0),size(fit.specs,2));
par = fit.param;
if isfield(par,'BL')
    par = rmfield(par,'BL');
end
doBaseline = false;
if isfield(par,'baseline')
    doBaseline = true;
    baseline = par.baseline;
end
isoH = 0;
if nargout > 2
    Pisos = [];
    Hisos = [];
    parent = [];
end
% Reverse direction of i??
for i = 1:5
    
    ind = round(fit.mz) == M-i;
    iM = par.massRange == M-i;
    if ~any(ind)
        continue;
    end
    
    if approx
        % If asked to approximate the impact of isotopes on the previous
        % mass, the function calls itself, without approximating for each of the masses before.
        % Not a great approximation, only intended for plotting.
        Yiso = getIsoYs(fit,compList,M-i,false);
        fit.specs(ind,:) = fit.specs(ind,:) - Yiso;
    end
    
    x = fit.mz(ind);
    iPeaks = round(peakList) == M-i;
    par.peakList = peakList(iPeaks);
    if isfield(fit,'deltaMassCal')
        par.dmfit = fit.deltaMassCal(iM);
    end
    if doBaseline
        par.baseline = baseline(ind,:);
    end
    % if isfield(par,'BL')
    %     find(par.massRange == M-i,1)
    %     par.BL
    %     par.BL = par.BL(par.massRange == M-i);
    % end
    
    if ~(sum(iPeaks) == 0)
        nams = compNams(iPeaks);
        if useH
            H = Hgen(iPeaks,:);
            H(isnan(H(:,1)),:) = [];
        else
            y = fit.specs(ind,:);
            [~,~,~,H,~,~] = fitPeaks(x, y, [], par);
        end
        Hiso = zeros(size(H,1)*10,size(H,2));
        Piso = nan(size(H,1)*10,1);
        if nargout > 2
            pnt = strings(size(H,1)*10,1);
        end
        jj = 1;
        for j = 1:sum(iPeaks)
            switch char(nams(j))
                case 'unknown'
                    continue
                otherwise
                    if ~isnan(str2double(nams(j)))
                        continue
                    end
                    %newIsotopes = tof_exact_isotope_masses(char(nams(j)),10);
                    newIsotopes = getIsoMass(nams(j),1e-7);
                    indH = round(newIsotopes(:,1))==M;
                    nIso = sum(indH);
                    if nIso > 0
                        Piso(jj:jj+nIso-1) = newIsotopes(indH,1);
                        Hiso(jj:jj+nIso-1,:) = H(j,:).*(newIsotopes(indH,2)./newIsotopes(1,2));
                        if nargout > 2
                            pnt(jj:jj+nIso-1) = nams(j);
                        end
                        jj = jj + nIso;
                    end
            end
        end
        rm = isnan(Piso) | sum(Hiso,2) == 0;
        Piso(rm) = [];
        Hiso(rm,:) = [];
        isoH = isoH + sum(Hiso,'all','omitnan');

        if nargout > 2
            pnt(rm) = [];
            Pisos = [Pisos; Piso];
            Hisos = [Hisos; Hiso];
            parent = [parent; pnt];
        end
        
        if isoH > 0
            Y = mkFigYs(x0, Piso, Hiso, par);
            isoYs = isoYs + squeeze(sum(Y,2));
%             figure;
%             plot(x0,squeeze(sum(Y,2)))
%             title(string(M-i))
        end
    end
end

end

