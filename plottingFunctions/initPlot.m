function [f,ax] = initPlot(sz,varargin)
%INITPLOT Initialize plot according to personal preference

if nargin == 0
    f = mkFigure();
    ax = mkAxes(f);
elseif isa(sz,"double")
    f = mkFigure(sz,varargin);
    ax = mkAxes(f,false);
else
    ax = mkAxes(f);
    f = mkFigure([{sz},varargin{:}]);
end
end

