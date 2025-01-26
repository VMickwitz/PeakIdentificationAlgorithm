
% specified by elemental symbol as array [number, precise mass, abundance]

% Values from : https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
lines = readlines('isoTable.txt');
nlines = length(lines);
nam = "";
arr = [];
isoStruct = struct;

for i = 2:nlines
    line = lines(i);
    if contains(line,'_')
        continue;
    end
    parts = strsplit(line,' ');
    
    if length(parts) < 4
        % Skip if not natural isotope
        continue;
    end
    if ~strcmp(parts(1),"")
        if ~strcmp(nam,"") && ~isempty(arr)
            isoStruct.(nam) = arr;
            arr = [];
        end
        parts = parts(2:end);
        nam = parts(1);
    elseif strcmp(nam,"H")
        parts = parts(2:end);
    end
    if strcmp(parts(4),"")
        % Skip if not natural isotope
        continue;
    end
    num = str2double(parts(2));
    mass = str2double(regexprep(parts(3),"\(\d*\)",""));
    abundance = str2double(regexprep(parts(4),"\(\d*\)",""));
    if isnan(abundance) || abundance > 1
        continue;
    end
    arr = [arr ;[num, mass, abundance]]; %#ok<AGROW>
    % if i > 20
    %     break
    % end
end

% Explicitly save values for certain number of common elements.
commonEl = ["H" "C" "N" "O" "F" "S" "I" "Br"];
ul = [60 40 20 30 30 20 20 10];

for i = 1:length(commonEl)
    el = commonEl(i);
    for j = 2:ul(i)
        nam = sprintf("%s%i",el,j);
        isoStruct.(nam) = getIsoMass(nam);
    end
end

isoStruct
isoStruct.H
isoStruct.C
isoStruct.Fe
isoStruct.U
isoStruct.C20
save("isoStruct","isoStruct")

%%
% clear;
% clearvars -global;

tic
iso1 = getIsoMass("N")
toc

tic
iso2 = tof_exact_isotope_masses('N')
toc

l = min(length(iso1),length(iso2));
t1 = sum(iso1(1:l,2))-1;
t2 = sum(iso2(1:l,2))-1;
dt = t1-t2
deltaP = iso1(1:l,2)-iso2(1:l,2);
plot(deltaP)
%%
fprintf("___________________________________________\n")
for i = 1:min(10,l)
    fprintf("m = %9.5f\n",iso2(i,1)-0.000548579909)
    fprintf("I1: %8.7f\n",iso1(i,2)./iso1(1,2)*100)
    fprintf("I2: %8.7f\n\n",iso2(i,2)./iso2(1,2)*100)
end