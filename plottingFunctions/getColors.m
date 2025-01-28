function [colors] = getColors(name,n)
%GETCOLORS Returns different colormatrices.
%   n is an additional input for grayscale.

blue = [0.14, 0.14, 0.8];
red = [0.86, 0.14, 0.14];
green = [0.27, 0.6, 0.21];
yellow = [1 0.83 0.05];
purple = [0.65, 0.3, 0.78];
lblue = [0.3, 0.75, 0.94];
winered = [0.68, 0.06, 0.22];
orange = [0.95, 0.54, 0.05];
lgreen = [0.52, 0.92, 0.45];
magenta = [0.92, 0.45, 0.82];

switch name
    case "default"
        colors = [0, 0.4470, 0.7410;
        0.8500, 0.3250, 0.0980;
        0.9290, 0.6940, 0.1250;
        0.4940, 0.1840, 0.5560;
        0.4660, 0.6740, 0.1880;
        0.3010, 0.7450, 0.9330;
        0.6350, 0.0780, 0.1840];
    case "mycolors"
        colors = [blue;red;green;yellow;purple;lblue;winered;lgreen;orange;magenta];
    case "poster"
        colors = [0, 0, 0;
            252, 163, 17;
            140, 140, 140]/255;
    case "grays"
        if nargin > 1
            colors = grayscale(n);
        else
            colors = grayscale(5);
        end
    otherwise
        colors = [0, 0.4470, 0.7410;
        0.8500, 0.3250, 0.0980;
        0.9290, 0.6940, 0.1250;
        0.4940, 0.1840, 0.5560;
        0.4660, 0.6740, 0.1880;
        0.3010, 0.7450, 0.9330;
        0.6350, 0.0780, 0.1840];
end

