function [idx] = checkMass(numel,clusters,seeds)
%CHECKMASS Checks if the compound can be made up of the clusters specified.
idx = false([size(numel,1),1]);
idel = false(size(seeds,1));

% Check if any molecule is a seed.
for i = 1:size(seeds,1)
    if any(all(numel<seeds(i,:),1),2) 
        % Remove obsolete seeds.
        idel(i) = true;
    else
        idx(all(numel==seeds(i,:),2)) = true;
    end
end
seeds(idel,:) = [];

idel = false(size(clusters,1));
for i = 1:size(clusters,1)
    % Remove obsolete clusters
    if any(all(numel<clusters(i,:),1),2)
        idel(i) = true;
    end
end
clusters(idel,:) = [];

for i = 1:size(clusters,1)
    % Remove cluster if possible recheck.
    ind = ~idx & all(numel >= clusters(i,:),2);
    idx(ind) = checkMass(numel(ind,:)-clusters(i,:),clusters,seeds);
end

end
