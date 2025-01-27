function [strOut] = unformat(str)
%UNFORMAT Summary of this function goes here
%   Detailed explanation goes here
str = replace(str,'\','\\\');
str = replace(str,'_','\_');
strOut = replace(str,'^','\^');
end

