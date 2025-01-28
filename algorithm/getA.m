function fit = getA(fit,doPeaks,doPlot)
% GETA Determines a suitable value of A.
% Method improved to be more robust than the one described in 
% Mickwitz et al. (2025).

Afit = fit.fitPrelim;
param = fit.param;

if nargin < 2 || isempty(doPeaks)
    doPeaks = false;
end

if nargin < 3 || isempty(doPlot)
    doPlot = true;
end

if isfield(param,'avgPeakLims')
    plims = param.avgPeakLims*length(param.massRange);
else
    plims = param.peakRange*length(param.massRange);
end

if isfield(Afit,'peakN') && ~doPeaks
    nPeaks = Afit.nPeaks;
    As = Afit.Avec;
    chis = Afit.chis;
    chiend = Afit.chiend;
else
    fprintf("Getting peak numbers as a function of A...\n")
    k = size(Afit.chi,2);
    inan = ~isnan(Afit.chi);
    chiend = nan(k,1);
    for i = 1:k
        %ind = find(inan(:,i);
        chiend(i) = Afit.chi(find(inan(:,i),1,'last'),i);
    end

    % Get all lower limits for A intervals
    %lims = mappingFunction(permute(Afit.Alims,[3,1,2]),fit.param.W);
    lims = permute(Afit.Alims,[3,1,2]).*chiend;
    lims = permute(lims,[2 3 1]);
    llim = squeeze(lims(:,1,:));
    llim = llim(:);
    % Get indices ordered from lowest to highest
    [Lsort,iSort] = sort(llim);
    % Remove nan and 0 values
    irm = isnan(Lsort) | Lsort==0;
    iSort(irm) = [];
    % Get the relevant values
    As = llim(iSort);
    % Add edge values for nicer plots, and initialization
    As = [As(1)*0.9; As];
    % Obtain indices of A values
    n = length(As);
    m = length(param.massRange);
    if doPeaks
        peakvec = repmat(squeeze(Afit.peaks(1,:,:)),[1,1,n]);
        H = squeeze(sum(Afit.H,2,'omitnan'));
        Hvec = repmat(squeeze(H(:,1,:)),[1,1,n]);
    end
    % Prepare a matrix for getting the second dimension index:
    massInd = repmat(1:m,[size(Afit.chi,1),1]);
    % Initialize
    nPeaks = repmat(min(fit.param.peakRange),[n,m]);
    chis = repmat(Afit.chi(1,:),[n,1]);
    % For every value of A, update the value at the updating mass only.
    for i = 2:n
        k = iSort(i-1); % Index for i:th value of As
        j = massInd(k); % 2nd dimension index
        % Update for this, and all coming values of As
        nPeaks(i:n,j) = Afit.numbers(k);
        chis(i:n,j) = Afit.chi(k);
        if doPeaks
            % Collect all peaks and heights for each value in As.
            ind = rem(k-1,13)+1;
            peakvec(:,j,i:n) = repmat(Afit.peaks(ind,:,j)',[1,1,n-i+1]);
            Hvec(:,j,i:n) = repmat(H(:,ind,j),[1,1,n-i+1]);
        end
    end
    % Save results:
    fit.fitPrelim.Avec = As;
    fit.fitPrelim.nPeaks = nPeaks;
    fit.fitPrelim.chis = chis;
    fit.fitPrelim.chiend = chiend;
    if doPeaks
        fit.fitPrelim.peakvec = peakvec;
        fit.fitPrelim.Hvec = Hvec;
    end
end

%chimu = mean(chis,2);
%chimu = mean(chis,2);
chimat = chis./chiend';
chimu = mean(chimat,2);
peakN = sum(nPeaks,2);

% Rvec = corr(chimat',param.W);
% Rpm = cumsum((Rvec>0)-(Rvec<0));
% Rs = smooth(Rvec,201);
% 
% start = find(peakN>(peakN(end-1)-4*m),1,'first');
ill = find(peakN > plims(1),1,'first');
iul = find(peakN < plims(2),1,'last');
init = [1, 2.2];

% Less than 2% of masses should reach the cap on peak number.
isMax = sum(nPeaks==nPeaks(end,:),2)/m;
iCap = find(isMax < 0.02,1,'last');
iul = min(iul,iCap);

try
    [iA, par] = getiA(peakN,chimu,init(1),m,ill,iul);
    [~,fval] = fitf(par,peakN,chimu);
catch me
    % This section is prone to errors. Try to find all the bugs.
    [~,ax] = initPlot();
    plot(As,peakN/m)
    ax.XScale = 'log';
    title("Error identifying A")
    xlabel("Parameter A")
    ylabel("Number of peaks per unit mass")
    rethrow(me)
end

if isempty(iA)
    error("iA empty. Check indexing.")
end
% Invert mapping of A (turns A into vector of length(W))
%fit.param.A = mappingFunction(As(iA),fit.param.W,true);
fit.param.A = As(iA)./chiend;

fit.fitPrelim.Aprime = As(iA);

if doPlot

    t = mkTiles([2 2]);
    ax = addTile(t);
    grid minor
    plot(As,peakN/m)
    ax.XScale = 'log';
    xlim([As(1), As(end)])
    ylabel("Number of peaks per unit mass")
    xlabel("Parameter A")
    hold on
    yyaxis right
    plot(As,isMax)
    hold off
    xline(As(iA),'k-','LineWidth',1.5);
    legend(["Average number of peaks" "Fraction at max"],'Location','northwest')
    
    iplot = peakN <= par(2);
    ax = addTile(t);
    grid minor
    plot(peakN/m,chimu);%./(peakN/peakN(end-1)))
    plot(peakN(iplot)/m,fval(iplot),'-.')
    xline(peakN(iA)/length(fit.param.massRange),'k-','LineWidth',1.5)
    ylabel("Mean relative $\chi^2$ value")
    xlabel("Number of peaks per unit mass")
    ax.YScale = 'log';
    ylim([1, 1+5*(chimu(iA)-1)])
    xlim([0,12])
    legend(["Data" "Fit"])

    ax3 = addTile(t);
    [AA,MM] = meshgrid(As,fit.param.massRange);
    colormap("parula");
    surface(ax3,AA,MM,nPeaks','EdgeColor','none');
    ax3.XScale = 'log';
    xlabel("Parameter A")
    ylabel("m/z")
    xline(As(iA),'k-')
    xlim([As(1), As(end)])
    cb = colorbar(ax3,'eastoutside');
    cb.Label.String = "Number of peaks";
    completeTiles();

    % addLabels(t,'northwest','alpha');
end

end

function [score, yhat] = fitf(par,x,y)
    try
    yhat = par(1)*(par(2)-x).^2;
    catch ME
        % length(x)
        % length(y)
        % any(isnan(x))
        % any(isnan(y))
        % par
        rethrow(ME)
    end
    score = mean((log(y)-yhat).^2);
    if nargout > 1
        yhat = exp(yhat);
    end
end

function [iA,par] = getiA(x,y,init,m,ill,iul)

iend = find(log(y)<1e-6,1,'first');
if isempty(iend)
    iend = length(x);
end
x = x(1:iend);
y = y(1:iend);
iul = min(iul,iend);
stepsz = 5;
iAs = ill:stepsz:iul;
n = length(iAs);
scores = nan(n,1);
bs = nan(n,2);

ll = [0, (x(end)-x(ill))/2];
ul = [init(1)*1000, x(end)+m/2];
init = [init, x(end)+m/8];

for i = 1:n
    iA = iAs(i);
    [par,scores(i)] = boundfmin(@(b)fitf(b,x(iA:end),y(iA:end)),init,ll,ul);
    bs(i,:) = par;
    init = par;

end

delta = bs(:,1).*bs(:,2).^2 - log(y(1));
iuse = delta < 0;

ind = find(diff(scores(iuse))>=0,1,'last');

% initPlot();
% plot(diff(scores(iuse)))

if isempty(ind)
    ind = find(diff(scores)>=-1e-6,1,'first');
    if isempty(ind)
        warning("No local minimum found in fitScore. Using highest slope value.")
        nvec = 1:n;
        inds = (nvec(iuse))*10;
        % x(inds)
        % scores(iuse)
        [~,ind] = max(diff(scores(iuse))./diff(x(inds)));
    end

end

ind = find(cumsum(iuse)==ind,1,'first');
par = bs(ind,:);

[~,yhat] = fitf(par,x,y);

% Needs improving
iA = find(log((y)./yhat)<5e-2,1,'first');
if iA > iul
    iA = iul;
    warning("A value likely too high.")
elseif iA < ill
    iA = ill;
    warning("A value likely too low.")
end

end
