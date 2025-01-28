function h = addTile(t,span)
%ADDTILE Summary of this function goes here
%   Detailed explanation goes here

if nargin == 1
    nexttile;
else
    nexttile(span);
end

h = mkAxes(t,false);

end
