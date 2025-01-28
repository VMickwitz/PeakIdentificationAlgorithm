function [f] = mkFigure(sz,varargin)
%PLOT0 Makes standard plot according to personal preference.

fsz = 14;
tfszm = 1.5;
axfszm = 1.4;
%fpos = [680 558 580 500]; % Pixels
doMenuBar = false;
if doMenuBar
    top_margin = [0 78]*0.75;
else
    varargin = [varargin{:} {'MenuBar', 'none'}];
    top_margin = [0 30]*0.75;
end
bot_margin = [0 30]*0.75;
margin = top_margin+bot_margin;
screensize = get(0,"ScreenSize")*0.75-[1,1,margin];
screensize = screensize(3:4);
figsz = [410,360];
%fpos = [510 418 410 360];
lineWidth = 2;

%fpos = pos+[0 0 expIns(1) + expIns(3) expIns(2)+expIns(4)];
if nargin == 0
    fpos = [(screensize-figsz)./4, figsz];
    f = figure('Units','points','Position',fpos,'NextPlot','add','PaperUnits','points');
else
    if isa(sz,"string") || isa(sz,"char")
        varargin = [{sz}, varargin{:}];
    else
        figsz = sz;
    end
    if any(figsz>screensize)
        figsz(figsz>screensize) = screensize(figsz>screensize);
    end
    fpos = [bot_margin+(screensize-figsz)./4, figsz];
    %fpos = [0,23,figsz];
    if isempty(varargin)
        f = figure('Units','points','Position',fpos,'NextPlot','add','PaperUnits','points');
    else
        f = figure('Units','points','Position',fpos,'NextPlot','add','PaperUnits','points',varargin{:});
    end
    %f = figure('Units','points','Position',fpos,'NextPlot','add','PaperUnits','points',varargin{:});
end
set(f,'DefaultLineLineWidth', lineWidth);
set(f,'DefaultAxesFontSize',fsz);
set(f,'DefaultAxesTickLabelInterpreter','latex');
set(f,'DefaultLegendInterpreter','latex');
set(f,'DefaultTextInterpreter','latex');
set(f,'DefaultAxesLabelFontSizeMultiplier',axfszm);
set(f,'DefaultAxesTitleFontSizeMultiplier',tfszm);
set(f,'DefaultLegendFontSize',fsz);
set(f,'DefaultLegendFontSizeMode','manual');
set(f,'DefaultLegendLocation','best');
set(f,'DefaultTextFontSize',fsz);
set(f,'DefaultAxesFontName','Times New Roman');
end

