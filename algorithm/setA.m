function [fit] = setA(fit,A,doPlot)
%SETA Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    doPlot = false;
end

fit.param.A = A./fit.fitPrelim.chiend;
fit.fitPrelim.Aprime = A;

if doPlot

    As = fit.fitPrelim.Avec;
    nPeaks = fit.fitPrelim.nPeak;
    chiend = fit.fitPrelim.chiend;
    chis = fit.fitPrelim.chis;
    chimat = chis./chiend';
    chimu = mean(chimat,2);
    m = length(fit.param.massRange);
    peakN = sum(nPeaks,2);
    isMax = sum(nPeaks==nPeaks(end,:),2)/m;
    iA = find(As < A,1,'last');

    t = mkTiles([2 1]);
    %f = t.Parent;
    %f.Position = [0 0 900 300];
    ax = addTile(t);
    grid minor
    plot(As,peakN/m)
    ax.XScale = 'log';
    xlim([1e-6 1e8])
    ylabel("Number of peaks per unit mass")
    xlabel("Parameter A")
    hold on
    yyaxis right
    plot(As,isMax)
    hold off
    xline(A,'k-','LineWidth',1.5);
    legend(["Average number of peaks" "Fraction at max"])

    ax = addTile(t);
    grid minor
    plot(peakN/m,chimu);
    plot(peakN/m,fval,'-.')
    xline(peakN(iA)/length(fit.param.massRange),'k-','LineWidth',1.5)
    ylabel("Mean relative $\chi^2$ value")
    xlabel("Number of peaks per unit mass")
    ax.YScale = 'log';
    xlim([0,12])
    legend(["Data" "Fit"])
    completeTiles();
    addLabels(t,'northwest','alpha');

end

end

