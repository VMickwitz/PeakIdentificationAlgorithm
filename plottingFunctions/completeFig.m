function completeFig(fh,ax)
%COMPLETEFIG Summary of this function goes here
%   Detailed explanation goes here

if nargin == 0
    fh = gcf;
    ax = gca;
elseif nargin == 1
    ax = fh.Children(1);
end
%fullscreen(fh,ax)
margin = [0.02 0.02];
margin = margin.*ax.Position([2 3]);
ax.Position([1 2]) = ax.Position([1 2]) - ax.OuterPosition([1 2])+margin;
fh.Position([3,4]) = ax.OuterPosition([3 4])+margin*2;
end

