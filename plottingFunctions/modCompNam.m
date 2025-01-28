function [latexComp] = modCompNam(comp)
%MODCOMPNAM Modifies compound name to latex compatible.
%   C10H12O14F0- becomes C$_{10}$H$_{12}$O$_{14}$$^{-}$

comp = string(comp);
if length(comp) > 1
    latexComp = strings(size(comp));
    for i = 1:length(comp)
        latexComp(i) = modify(comp(i));
    end
else
    latexComp = modify(comp);
end

end

function name = modify(formula)
if strcmp(formula,"unknown")
    name = "unknown";
elseif ~isnan(str2double(formula))
    name = formula;
else
    name = regexprep(formula,"\s","");
    name = regexprep(name,"\D+[0]","");
    name = regexprep(name,"\d*","$_{$&}$");
    name = regexprep(name,"[+-]\d*","$^{$&}$");
    order = ["C" "H" "O" "N" "F" "S" "Cl"];
    [parts,el] = regexp(name,"([a-zA-Z]+)[$]_{\d*}[$]|[$][^]{[+-]\d*}[$]","match","tokens");
    idx = nan(1,size(parts,2));
    nel = length(el);
    range = 1:length(order);
    for i = 1:nel
        if isempty(el{i})
            continue
        end
        try
            ind = strcmp(order,el{i});
        catch ME
            order
            el{i}
            el
            rethrow(ME)
        end
        if any(ind)
            idx(i) = range(ind);
        end
    end
    [~,isort] = sort(idx,2,"ascend","MissingPlacement","last");
    parts = parts(isort);
    name = strjoin(parts,"");
    name = regexprep(name,"[$][$]","");
end
end