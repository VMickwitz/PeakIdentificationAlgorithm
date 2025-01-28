function [t] = mkTiles(arrangement,varargin)
%MKTILES Summary of this function goes here
%   Detailed explanation goes here

if nargin == 0 || isempty(arrangement)
    args = [{"flow",'Padding','tight','TileSpacing','compact'} varargin{:}];
elseif length(arrangement) == 2
    args = [{arrangement(1),arrangement(2),'Padding','tight','TileSpacing','compact'} varargin{:}];
else
    args = [{arrangement,'Padding','tight','TileSpacing','compact'} varargin{:}];
end

%fpos = pos+[0 0 expIns(1) + expIns(3) expIns(2)+expIns(4)];

f = mkFigure;
%f.Position = [0 0 1385 720];
f.Position = [0 0 1280 642];
f.Visible = 'on';
fsz = get(f,'DefaultAxesFontSize');
axfszm = get(f,'DefaultAxesLabelFontSizeMultiplier');
interpreter = get(f,'DefaultTextInterpreter');
tfszm = get(f,'DefaultAxesTitleFontSizeMultiplier');
lineWidth = get(f,'DefaultLineLineWidth');

%t = tiledlayout(f,arrangement{:},'Padding','tight','TileSpacing','compact',varargin{:});
t = tiledlayout(f,args{:});

t.XLabel.Interpreter = interpreter;
t.YLabel.Interpreter = interpreter;
t.Title.Interpreter = interpreter;
t.PositionConstraint = "innerposition";
%t.PositionConstraint = "outerposition";
t.XLabel.FontSize = fsz*axfszm;
t.YLabel.FontSize = fsz*axfszm;
t.Title.FontSize = fsz*2;
set(t,'DefaultLineLineWidth', lineWidth);
set(t,'DefaultAxesTickLabelInterpreter','latex');
set(t,'DefaultAxesLabelFontSizeMultiplier',axfszm);
set(t,'DefaultLegendInterpreter','latex');
set(t,'DefaultTextInterpreter','latex');
set(t,'DefaultAxesLabelFontSize',axfszm);
set(t,'DefaultAxesTitleFontSizeMultiplier',tfszm);
set(t,'DefaultLegendFontSize',fsz);
set(t,'DefaultLegendFontSizeMode','manual');
set(t,'DefaultLegendLocation','best');
set(t,'DefaultTextFontSize',fsz);
end

