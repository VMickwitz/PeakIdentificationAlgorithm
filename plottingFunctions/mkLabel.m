function [labelStr] = mkLabel(name,unit)
%MKLABEL Summary of this function goes here
%   Detailed explanation goes here

unit = fmtUnit(unit);
labelStr = sprintf("%s (%s)",name,unit);
end

function unit = fmtUnit(ustr)
% Changes barebones unit notation to latex compatible.
    unit = regexprep(ustr,"-?\d+","$^{$&}$");
end