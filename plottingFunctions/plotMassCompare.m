function [tl] = plotMassCompare(fit,M,list,separate)
%PLOTMASSCOMPARE Summary of this function goes here
%   Detailed explanation goes here

if nargin < 4
    separate = false;
end

if nargin > 2
    if ~separate
        tl = mkTiles([2 2]);
    end
    doManual = true;
    Rlegend = ["Assigned" "Free fit" "Manual"];
    %Rlegend = ["Assigned" "Free fit" "Generated"];
    rep = [3 1];
else
    if ~separate
        tl = mkTiles([1 3]);
    end
    doManual = false;
    Rlegend = ["Assigned" "Free fit"];
    rep = [2 1];
end

yls = [0 0];
prefit = getPrelimFit(fit,M);
% Make BL the same for all figures:
BL = prefit.param.BL;
fit.param.BL(:) = BL;

if separate
    [~,a] = initPlot();
else
    a = addTile(tl,1);
end
plotMass(fit,M,'axes',a,'label',Rlegend(1));
%set(gca,'YScale','log')

peaks = fit.peaks(:,fit.param.massRange==M);
if ~all(isnan(peaks))
    xls = [min(peaks)-0.15, max(peaks)+0.1];
else
    xls = [-inf, inf];
end
yls = updateLimits(yls);

yyaxis right;
res = gca().Children(2).YData;
x = gca().Children(2).XData;
res = repmat(res,rep);

if separate
    completeFig();
    [~,b] = initPlot();
else
    b = addTile(tl,2);
end
plotMass(prefit,M,'axes',b,'label',Rlegend(2));
yls = updateLimits(yls);

yyaxis right;
res(2,:) = gca().Children(2).YData;
%rls = updateLimits(rls);

if doManual
    if separate
        completeFig();
        [~,c] = initPlot();
    else
        c = addTile(tl,3);
    end
    manualFit = getListFit(fit,list,M);
    manualFit.param.BL(:) = BL;
    plotMass(manualFit,M,'axes',c,'label',Rlegend(3));
    yls = updateLimits(yls);

    yyaxis right;
    res(3,:) = gca().Children(2).YData;
    %rls = updateLimits(rls);
    if separate
        completeFig();
        [~,d] = initPlot();
    else
        d = addTile(tl,4);
    end
else
    if separate
        completeFig();
        [~,d] = initPlot();
    else
        d = addTile(tl,3);
    end
end
yls(1) = -0.02*yls(2);
a.YAxis(1).Limits = yls;
b.YAxis(1).Limits = yls;
if doManual
    c.YAxis(1).Limits = yls;
end
grid minor
plot(d,x,res)
xlabel("m/z (Th)")
ylabel("Residual")
legend(Rlegend,'Location','northwest')
title('Residuals','Units','normalized','Position',[0.83,0.85,0])
rls = gca().YLim;

a.XLim = xls;
a.YLim = rls;
b.XLim = xls;
b.YLim = rls;
if doManual
    c.XLim = xls;
    c.YLim = rls;
end
% I have no idea why this works....


% Align grids on axes.
la = a.YAxis(1);
ra = a.YAxis(2);
lticks = la.TickValues;
rticks = ra.TickValues;
nticks = numel(lticks);
if nticks ~= numel(rticks)
    if nticks > numel(rticks)
        delta = rticks(2)-rticks(1);
    else
        delta = (rticks(2)-rticks(1))*2;
    end
    if rem(nticks,2) == 0
        min(res,[],'all')
        max(res,[],'all')
        if abs(min(res,[],'all')) > max(res,[],'all')
            rtickNew = -(delta*(nticks/2)):delta:(delta*(nticks/2-1));
        else
            rtickNew = -(delta*(nticks/2-1)):delta:(delta*(nticks/2));
        end
    else
        rtickNew = -(delta*((nticks-1)/2)):delta:(delta*((nticks-1)/2));
    end
else
    rtickNew = rticks;
end

d1 = lticks(end)-lticks(1);
d2 = rtickNew(end)-rtickNew(1);
rls(1) = rtickNew(1)-d2*(lticks(1)-yls(1))/d1;
rls(2) = rtickNew(end)+d2*(yls(2)-lticks(end))/d1;

a.YLim = rls;
a.YAxis(2).TickValues = rtickNew;
b.YLim = rls;
b.YAxis(2).TickValues = rtickNew;
if doManual
    c.YLim = rls;
    c.YAxis(2).TickValues = rtickNew;
end
d.YLim = rls;
d.YTick = rtickNew;



if separate
    completeFig();
else
    completeTiles();
end

end

function [ls] = updateLimits(ls)
ils = gca().YLim;
if ils(1)<ls(1)
    ls(1) = ils(1);
end
if ils(2)>ls(2)
    ls(2) = ils(2);
end
end

