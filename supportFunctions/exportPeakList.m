function [success] = exportPeakList(fit,filepath,fmt,varargin)
%EXPORTPEAKLIST Summary of this function goes here
%   Detailed explanation goes here

success = false;
nvararg = length(varargin);

if nargin < 2 || isempty(filepath)
    [fname,fpath] = uiputfile("*.txt");
    filepath = [fpath fname];
else
    fname = regexp(filepath,"\w+[.]txt$",'match');
end

if nargin < 3 || isempty(fmt)
    fmt = "tofware";
end

switch fmt
    case "tofware"
        fmt = 1;
    case "toftools"
        fmt = 2;
        error("tofTools format not yet implemented.")
    otherwise
        error("Unknown export format.")
end

incUnknown = true;
fitUnknown = true;
defUnknown = true;
for i = 1:2:nvararg
    switch varargin{i}
        case 'includeUnknown'
            incUnknown = varargin{i+1};
        case 'fitUnknown'
            fitUnknown = varargin{i+1};
        case 'defUnknown'
            defUnknown = varargin{i+1};
        otherwise
            warning("Unknown argument for exportPeakList: %s",varargin{i})
    end
end


fID = fopen(filepath,'wt');
listlen = length(fit.PeakList);

h = waitbar(0,'Formatting data for export...');
masses = fit.PeakList;
names = erase(fit.CompNams," ");
charges = regexp(names,'[+-]\d*','match');
nams = erase(modNams(fit.CompNams,[],false)," ");

waitbar(0,h,sprintf("Writing data to %s...",unformat(fname)));
unknownCounter = 1;
if fmt == 1
    fprintf(fID,"def\tfit\tion\tx0\ttag\tsumFormula\tx_Lo\tx_Hi\td2_Ctr\td2_Lo\td2_Hi\td1_Ctr\td1_Lo\td1_Hi	d0_Ctr\td0_Lo\td0_Hi\tcalFac\tcalUnit\tcharge\tionizFrac\tfragOf\tisotopeOf\tcomments\tstandards\tfamilies\tbaseline\tintCalIn\n");

    for i = 1:listlen
        if isempty(charges{i})
            charge = "";
        elseif strcmp(charges{i},"-")
            charge = "-1";
        elseif strcmp(charges{i},"+")
            charge = "+1";
        else
            charge = "";
        end

        isUnknown = strcmp("unknown",names(i));
        if isUnknown
            names(i) = sprintf("<unknown%04i>",unknownCounter);
            unknownCounter = unknownCounter+1;
        end

        if incUnknown && defUnknown
            fprintf(fID,"%i\t%i\t%s\t%9.6f\t\t%s\t\t\t\t\t\t\t\t\t\t\t\t\t\t%s\t\t\t\t\t\t\t\t\n",1,fitUnknown,names(i),masses(i),nams(i),charge);
        else
            if isUnknown && incUnknown
                fprintf(fID,"%i\t%i\t%s\t%9.6f\t\t%s\t\t\t\t\t\t\t\t\t\t\t\t\t\t%i\t\t\t\t\t\t\t\t\n",~isUnknown,fitUnknown,names(i),masses(i),nams(i),charge);
            else
                fprintf(fID,"%i\t%i\t%s\t%9.6f\t\t%s\t\t\t\t\t\t\t\t\t\t\t\t\t\t%i\t\t\t\t\t\t\t\t\n",1,fitUnknown,names(i),masses(i),nams(i),charge);
            end
        end
        if rem(i,floor(listlen/100)) == 0
            waitbar(i/listlen,h);
        end
    end
end
fclose(fID);
delete(h);
success = true;
