function addLabels(t,location,type)
%ADDLABELS Adds labels to axes in tiled layout.

dist = 0.1; % [0 1]

if nargin < 2 || isempty(location)
    lmp = [1 0;0 1;0.5 0;0 1];
    vec = [0 -1];
else
    switch location
        case 'north'
            lmp = [1 0;0 1;0.5 0;0 1];
            vec = [0 -1];
        case 'south'
            lmp = [1 0;0 1;0.5 0;0 0];
            vec = [0 1];
        case 'east'
            lmp = [1 0;0 1;1 0;0 0.5];
            vec = [-1 0];
        case 'west'
            lmp = [1 0;0 1;0 0;0 0.5];
            vec = [1 0];
        case 'northeast'
            lmp = [1 0;0 1;1 0;0 1];
            vec = [-1 -1];
        case 'northwest'
            lmp = [1 0;0 1;0 0;0 1];
            vec = [1 -1];
        case 'southeast'
            lmp = [1 0;0 1;1 0;0 0];
            vec = [-1 1];
        case 'southwest'
            lmp = [1 0;0 1;0 0;0 0];
            vec = [1 1];
        otherwise
            error("Unknown Location")
    end
end


if nargin < 3 || isempty(type)
    type = 1;
else
    switch type
        case 'Alpha'
            type = 1;
        case 'alpha'
            type = 2;
        case 'num'
            type = 3;
        otherwise
            warning("Could not interpret type. Using alpha.")
            type = 2;
    end
end

if type == 1
    labels = char(65:90);
elseif type == 2
    labels = char(97:122);
elseif type == 3
    labels = string(1:28);
end

f = t.Parent;
scale = f.InnerPosition(3:4);
axs = t.Children;
unit = f.Units;
n = length(axs);
j = 1;
fsz = 18;
interpreter = 'latex';
p0 = [1 1];
p1 = [0 0];
isAx = false(n,1);
for i = 1:n
    isAx(i) = isequal(class(axs(i)),'matlab.graphics.axis.Axes');
    if isAx(i)
        if p0(1) > axs(i).InnerPosition(3)/scale(1)
            p0(1) = axs(i).InnerPosition(3)/scale(1);
            p1(1) = axs(i).InnerPosition(3);
        end
        if p0(2) > axs(i).InnerPosition(4)/scale(2)
            p0(2) = axs(i).InnerPosition(4)/scale(2);
            p1(2) = axs(i).InnerPosition(3);
        end
    end
end
p1(1:2) = min(p1(1:2));
pos = [0 0 0 0];

offSet = floor(dist.*p1);

for i = 1:n
    ind = n-i+1;

    if ~isAx(ind)
        continue;
    end
    % if ~isprop(axs(ind),'InnerPosition')
    %     continue;
    % end

    pos(1:2) = axs(ind).InnerPosition*lmp + (vec-1).*offSet;
    pos(3:4) = 2*offSet;
    pscale = pos./[scale scale];

    str = sprintf("(%s)",labels(j));
    tb = annotation(f,'textbox',pscale,'String',str,'FitBoxToText','off','EdgeColor','none','FontSize',fsz,'Interpreter',interpreter,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle','Units',unit);
    tb.Position = pos;
    j = j+1;

    % Show box size and corner vectors:
    % p0 = p1./scale;
    % pos2 = [axs(ind).InnerPosition*lmp./scale pos(3:4)]; % corners
    % arr = annotation(f,'arrow',[pos2(1) pscale(1)+dist*p0(1)],[pos2(2) pscale(2)+dist*p0(2)]);
    % %pos.*[scale scale]
    % rec = annotation(f,'rectangle',pscale,'Units',unit);
    % rec.Position = pos;
end

end

