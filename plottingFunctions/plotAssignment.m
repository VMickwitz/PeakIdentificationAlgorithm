function plotAssignment(mz,specs,peakList,compList,Hpeaks,Hcomps)
%ASSIGNMENTPLOT In progress plot of assignment process.


mz = fit.mz;
param = fit.param;
specs = fit.specs;
BL 
dmcal = 0;
y1s = mkFigYs(mz, [peakList,compList], [Hpeaks; Hcomp], param);
initPlot();
plot(mz,sum(specs,2),'k-')
plot(mz-dmcal,sum(y1s,2)+BL(ind),'rx')
plot(mz-dmcal,y1s)
for i = 1:length(param.peakList)
    xline(param.peakList(i),'k-','LineWidth',1.5)
end
for j = 1:length(compList)
    if j == 1
        xline(compList(j),'r--','LineWidth',1.5)
    else
        xline(compList(j),'k:','LineWidth',1.5)
    end
end

end
