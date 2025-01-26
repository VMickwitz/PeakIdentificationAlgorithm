load("isoStruct.mat")
nams = fieldnames(isoStruct);
massStruct = struct;

ind = regexp(nams,"\d+",'once');
for i = 1:length(nams)
    if isempty(ind{i})
        nams{i}
        massDat = isoStruct.(nams{i});
        [~,imax] = max(massDat(:,3));
        massStruct.(nams{i}) = massDat(imax,2);
    end
end
save("massStruct.mat","massStruct")