function [f] = completeTiles(t,margin,ax0)
%COMPLETETILES Run after adding all tiles to tiled layout.

if nargin == 0
    ax0 = get(gca);
    t = ax0.Parent;
    sz0 = ax0.InnerPosition(3:4);
    margin = [0 0];
elseif nargin < 2 || isempty(margin)
    sz0 = findMinSz(t);
    margin = [0 0];
else
    sz0 = ax0.InnerPosition(3:4);
end

sz = [340 240];
%sz = [272 200];
screenSz = [1385 720]-[0 5];
%screenSz = [960 460];

f = t.Parent;

t.Units = 'points';
f.Units = 'points';
t.OuterPosition(1:2) = 0;
deltaPos = t.OuterPosition(3:4)-t.InnerPosition(3:4);
if any(t.OuterPosition(3:4) > f.Position(3:4))
    t.Position(3:4) = f.Position(3:4);
end

newsz = -1;
while any(newsz < 0)
    gsz = t.GridSize([2 1]);
    pos0 = t.InnerPosition(3:4);
    padding = pos0 - gsz.*sz0;
    newsz = gsz.*sz+padding;
    if any((newsz+deltaPos+margin) > screenSz)
        % If figure larger than screen (any dimension), rescale)
        ratio = min(screenSz./(newsz+deltaPos+margin));
        newsz = newsz*ratio;
    end
    if any(newsz < 0)
        sz = sz.*0.8;
    end
end

t.InnerPosition(3:4) = newsz;
f.Position(3:4) = newsz+deltaPos+margin;
t.OuterPosition = [0 0 f.Position(3:4)]-[0 0 10 20];

drawnow;
end

function sz0 = findMinSz(t)

axs = t.Children;
n = length(axs);
sz0 = -1;
for i = 1:n
    if ~isprop(axs(i),'InnerPosition')
        continue;
    elseif sz0 == -1
        sz0 = axs(i).InnerPosition([3 4]);
    end
    if axs(i).InnerPosition(3) < sz0(1)
        sz0(1) = axs(i).InnerPosition(3);
    end
    if axs(i).InnerPosition(4) < sz0(1)
        sz0(2) = axs(i).InnerPosition(4);
    end
end

end