function [elements,numbers,charge] = getElements(names)
%Returns list of elements, numbers of atoms and charges.
%   (Make smaller faster version for single compound)


nNams = length(names);
names = cellstr(names);
charge = names{1};
charge = charge(end);
sz = 16;
nel = 0;
elements = string(nan(sz,1));
numbers = zeros(sz,nNams);

for i = 1:nNams
    nam = names{i};
    indchar = isletter(nam);
    modnam = nam;
    modnam2 = nam(1:end-1);
    modnam(~indchar) = ' ';
    modnam2(indchar) = ' ';
    while true
        [el,modnam] = strtok(modnam);
        [num,modnam2] = strtok(modnam2);
        if isempty(el)
            break
        end
        [elements,numbers,sz,nel] = addElement(elements,numbers,sz,nel,i,el,num);
    end
    
end

elements(nel+1:end) = [];
numbers(nel+1:end,:) = [];
end

function [elements,nums,sz,nel] = addElement(elements,nums,sz,nel,nn,el,num)
    ind = strcmp(elements,el);
    ok = ~any(ind);
    if ok
        nel = nel+1;
        ind = nel;
        if nel > sz
            elements = [elements; string(nan(sz,1))];
            nums = [nums; zeros(sz,nNams)];
            sz = sz*2;
        end
        elements(nel) = el;
    end
    nums(ind,nn) = str2double(num);
end