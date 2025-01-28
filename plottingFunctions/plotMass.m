function [f,a] = plotMass(fit,M,varargin)
%PLOTMASS Summary of this function goes here
%   Detailed explanation goes here


doLim = true;
legLoc = 'northwest';
label = "";
xls = [-inf, inf];
%yls = [-inf, inf];
doAx = false;
subLegend = false;
splitLegend = false;
dobf = false;
doLabels = true;
sparsel = false;
if nargin < 3
    varargin = {};
end
for i = 1:floor(length(varargin)/2)
    ind = 2*i-1;
    switch varargin{ind}
        case 'xls'
            xls = varargin{ind+1};
        case 'yls'
            yls = varargin{ind+1};
        case 'label'
            label = varargin{ind+1};
        case 'axes'
            ax = varargin{ind+1};
            doAx = true;
        case 'subLegend'
            subLegend = varargin{ind+1};
        case 'splitLegend'
            splitLegend = varargin{ind+1};
        case 'bf'
            dobf = true;
            bf = varargin{ind+1};
        case 'labels'
            doLabels = varargin{ind+1};
        case 'sparseLegend'
            sparsel = varargin{ind+1};
            if sparsel
                subLegend = varargin{ind+1};
            end
        otherwise
            warning("No argument called %s defined.",varargin{i})
    end
end

mz = fit.mz;
weights = fit.param.weights;
modnams = true;
if isfield(fit,'PeakList')
    iList = (round(fit.PeakList) == M);
    nams = fit.CompNams(iList);
    peaks = fit.PeakList(iList);
else
    peaks = fit.peaks;
    nams = strsplit(sprintf("%6.4f;",peaks),";");
    nams(end) = [];
    modnams = false;
end

if any(fit.peaks < xls(1))
    xls(1) = min(peaks)-0.05;
end
if any(fit.peaks > xls(2))
    xls(2) = max(peaks)+0.05;
end


atMass = round(mz) == M;
mz = mz(atMass);
ind = fit.param.massRange == M;

W = 1;%fit.param.W(ind);
specs = fit.specs(atMass,:) - fit.param.BL(ind).*weights*W;
%specs2 = fit.specsFinal(atMass,:);
H = fit.H(:,:,ind)*W;

H = H(~isnan(H(:,1)),:);

%iso1 = specs - specs2;
% list1.compNams = fit.CompNams;
% list1.peakList = fit.PeakList;
% isoFit = fit;
% isoFit.specs = fit.specsFinal;
% iso = getIsoYs(isoFit,list1,M);
iso = fit.isoSpecs(atMass,:);
plotIso = true;
if all(isnan(iso))
    plotIso = false;
end

[~, iSort] = sort(sum(H,2),'descend');
peaks = peaks(iSort);
nams = nams(iSort);
H = H(iSort,:);
if dobf
    bf = bf(iSort);
end

if isempty(peaks)
    y1s = nan(size(specs,1),1);
    yfit = zeros(size(specs));
else
    y1s = mkFigYs(mz, peaks, H, fit.param);
    yfit = squeeze(sum(y1s,2));% + fit.param.BL(fit.param.massRange == M).*weights;
end

pCol = getColors("poster");
xcol = pCol(2,:);
wdth = 2;

if modnams && ~isempty(nams)
    nam1 = modCompNam(nams');
    if any(strcmp(nam1,""))
        ifix = strcmp(nam1,"");
        nam1(ifix) = string(nams(ifix));
    end
else
    nam1 = nams;
end

if dobf
    for i = 1:length(nam1)
        if bf(i)
            nam1(i) = "\textbf{"+nam1(i)+"}";
        end
    end
end

if subLegend
    legStr1 = nam1;
elseif plotIso
    legStr1 = ["Signal" "Total fit" "Residual" "Isotopes" nam1];
else
    legStr1 = ["Signal" "Total fit" "Residual" nam1];
end

resColor = [0.83 0.83 0.83];

if doAx
    a = ax;
    f = ax.Parent;
else
    [f,a] = initPlot();
end
grid minor
%co = get(gca,'ColorOrder');
%set(gca,'ColorOrder',co([1 3 2 4 5 6 7],:)) % For poster
l0 = plot(a,mz,sum(specs,2),'k-','LineWidth',wdth);
hold on
l1 = plot(mz,sum(yfit,2)+sum(iso,2,'omitnan'),'x','Color',xcol,'LineWidth',wdth);
if plotIso
    l2 = plot(mz,sum(iso,2),'k:','LineWidth',wdth);
else
    set(gca,'ColorOrderIndex',4)
end
a.ColorOrderIndex = 1;
l3 = plot(mz,sum(y1s,3),'-','LineWidth',wdth);
if trapz(mz,sum(specs,2))<0
    %fprintf("BL failed?")
    %yline(fit.param.BL(fit.param.massRange == M).*sum(weights),'r:',...
    %'LineWidth',3);
    warning("BL mistake.")
end
hold off
yyaxis right
a.YAxis(2).Color = [0 0 0];
hold on
l4 = plot(mz,sum(specs,2)-sum(yfit,2)-sum(iso,2,'omitnan'),'-','LineWidth',1.5,'Color',resColor,'ZData',zeros(size(mz)));
plot(mz,zeros(size(mz)),'-','Color',resColor*0.9,'LineWidth',0.7,'ZData',zeros(size(mz))-1);
hold off
if doLabels
    ylabel("Residual")
end
yyaxis left
if doLabels
    ylabel("Signal")
end
nChil = length(a.Children);
for i = 1:nChil
    if isprop(a.Children(i),'Xdata')
        a.Children(nChil-i+1).ZData = i.*ones(size(a.Children(i).XData));
    else
        i
        a.Children(i)
    end
end
a.SortMethod = 'depth';
addPeaks(peaks, nams)
%set(gca,'FontSize',tickSz);
if subLegend
    if ~isempty(legStr1)
        if sparsel
            legend(l3,legStr1,'Location',legLoc,'Box',"off");
            %a.Legend.BackgroundAlpha = 0.1;
        else
            legend(l3,legStr1,'Location',legLoc);
        end
    end
elseif plotIso
    if splitLegend
        %a2 = axes('position',a.Position,'visible','off');
        legend(a,l3,legStr1(5:end),'Location','northwest')
        %legend(a,[l0; l1; l4; l2],legStr1(1:4),'Location','southwest')
    else
        legend([l0; l1; l4; l2; l3],legStr1,'Location',legLoc)
    end
else
    if splitLegend
        legend(a,l3,legStr1(4:end),'Location','northwest')
        %legend(a,[l0; l1; l4],legStr1(1:3),'Location','southwest')
        %a2 = axes('position',a.Position,'visible','off');
    else
        legend([l0; l1; l4; l3],legStr1,'Location',legLoc)
    end
end
if ~strcmp(label,"")
    title(label,'Units','normalized','Position',[0.83,0.85,0])
end

if doLabels
    xlabel("Mass to charge (Th)")
end
%ylabel("Signal")
%set(gca,'YTickLabel',[])
%a.XTickLabel(1) = [];
%set(gca,'XTick',xls(1)+0.05:0.1:xls(2)+0.05)
xlim(xls);
%ylim(yls);
%ylim([1e-5 inf])

%fullscreen();

end


function addPeaks(peaks,nams)

for i = 1:length(peaks)
    switch char(nams(i))
        case 'unknown'
            xline(peaks(i),'r--');
        otherwise
            xline(peaks(i),'k-');
    end
end
end
