function [success] = addAppPaths()
%ADDPATHS Adds algorithm files to matlab path
%   Detailed explanation goes here
success = true;
try
    cpath = mfilename("fullpath");
    if contains(cpath,"/")
        cpath = strrep(cpath,"/addAppPaths","");
    else
        cpath = strrep(cpath,"\addAppPaths","");
    end
    fpath = genpath(cpath);
    addpath(fpath);
catch
    success = false;
end
end
