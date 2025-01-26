function [modNams] = modNams(nams,order,incCharge)
% Changes compound names to match find_composition abbreviated names.
if nargin < 3
    incCharge = true;
end
if nargin == 1 || isempty(order)
    order = {'C' 'F' 'H' 'S' 'O' 'N' 'Br' 'Si' 'Zn' 'Cl' 'P'};
    %order = makeOrder(nams,order);
end

modify = @(nam) translate(nam,order,incCharge);
modNams = cellfun(modify,nams);

end

function newOrder = makeOrder(nams,order)
    elements = getElements(nams);
    inOrder = zeros(1,length(elements));
    for i = 1:length(order)
        inOrder(strcmp(order{i},elements)) = 1;
    end
    newOrder = horzcat(order,elements(~inOrder));
end

function [modNam, order] = translate(nam,order,incCharge)
if nam(end) == '+'
    charge = '+';
    nam = nam(1:end-1);
elseif nam(end) == '-'
    charge = '-';
    nam = nam(1:end-1);
else
    charge = [];
end

iChar = isstrprop(nam,'alpha');
iLC = ismember(nam,'a':'z');
indLC = find(iLC)-1; % Indices of elements with lower case letter
ind = find(iChar & ~iLC); % Indices of upper case letters
n = zeros(length(order),1);
l = sum(iChar)-sum(iLC);
for i = 1:l
    if any(ind(i)==indLC)
        cha = nam([ind(i) ind(i)+1]);
        if i == l
            if (ind(i)) == length(nam)
                num = '1';
            else
                num = nam(ind(i)+2:end);
            end
        else
            num = nam(ind(i)+2:ind(i+1)-1);
        end
    else
        cha = nam(ind(i));
        if i == l
            if (ind(i)) == length(nam)
                num = '1';
            else
                num = nam(ind(i)+1:end);
            end
        else
            num = nam(ind(i)+1:ind(i+1)-1);
        end
    end
    if isempty(regexp(num,'\d*','once'))
        num = '1';
    end
    spot = strcmp(cha, order);
    if ~any(spot)
        order = [order, cha]; %#ok<AGROW>
        n = [n;0]; %#ok<AGROW>
        spot = strcmp(cha, order);
    end
    n(spot) = n(spot) + str2double(deblank(num));
end
modNam = "";
for i = 1:size(order,2)
    if n(i) > 0
        modNam = sprintf("%s%s%i ",modNam,order{i},n(i));
    end
end
modNam = strtrim(modNam);
if ~isempty(charge) && incCharge
    modNam = sprintf("%s%c",modNam,charge);
end
if strcmp(modNam,"")
    testNam = str2double(nam);
    if isnan(testNam)
        modNam = "unknown";
    else
        modNam = string(nam);
    end
end

end