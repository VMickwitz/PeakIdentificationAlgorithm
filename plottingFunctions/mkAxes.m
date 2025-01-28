function [h] = mkAxes(f,usePos)
%MKAXES Summary of this function goes here
%   Detailed explanation goes here

pos = [0 0 340 240];
%pos = [0 0 272 220];
expIns = [60, 60, 10, 60];
if nargin < 2
    usePos = true;
end

h = newplot(f);
h.Units = 'points';
if usePos
    h.Position([1 2]) = expIns([1,2]);
    h.Position([3 4]) = pos([3 4]);
    h.PositionConstraint = "innerposition";
elseif ~strcmp(f.Type,'tiledlayout')
    h.Position([1 2]) = expIns([1 2]);
    h.Position([3 4]) = f.InnerPosition(3:4)-expIns(3:4)-expIns(1:2);
    h.PositionConstraint = "innerposition";
end
h.Box = 'off';
h.NextPlot = 'add';
h.LineWidth = 1;
%h.ColorOrder = getColors("mycolors");
h.ColorOrder = getColors("default");
grid on;
% grid minor;
axtoolbar(h,["zoomin", "zoomout", "pan", "restoreview", "export"]);
h.Toolbar.Visible = 'on';
end

