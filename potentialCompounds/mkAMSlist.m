clear;
mlim = [49.5 150.5];
el = ["C" "H" "O" "N" "S"];
charge = "+";
OtoC = [0 2];
HtoC = [0 4];
ul = [12 26 8 1 2];
ll = [0 0 0 0 0];

masses = [];
numel = [];
names = strings;

% Group1: Carbon, Hydrogen, Oxygen, Nitrogen
seeds = [1 0 0 0 0];
% Charged by cluster with nitrate or deprotonation (defluorization?)
ions = [0 0 0 0 0];
% single atoms
clusters = diag(ones(5,1));
nmax = [12 26 8 1 0];

seeds = addCharge(ions,seeds);

[mass,nums,name] = buildComp(el,seeds,clusters,nmax,mlim,charge);
masses = [masses; round(mass,4)];
numel = [numel; nums];
names = [names; name];

% Group2: with Sulfur
seeds = [0 0 0 0 1];
% Charged by cluster with nitrate or deprotonation (defluorization?)
ions = [0 0 0 0 0];
% single atoms
clusters = diag(ones(5,1));
nmax = [4 6 6 0 1];

seeds = addCharge(ions,seeds);

[mass,nums,name] = buildComp(el,seeds,clusters,nmax,mlim,charge);
masses = [masses; round(mass,4)];
numel = [numel; nums];
names = [names; name];

names = names(2:end);
OCratio = numel(:,3)./numel(:,1);
HCratio = numel(:,2)./numel(:,1);
C0 = numel(:,1) == 0;
irm = OCratio < OtoC(1) | OCratio > OtoC(2);
irm = irm | numel(:,2) < (numel(:,1)-4);
irm = irm | HCratio < HtoC(1) | HCratio > HtoC(2);
irm(C0) = false;
irm = irm | any(numel>ul,2);
irm = irm | any(numel<ll,2);

names(irm) = [];
numel(irm,:) = [];
masses(irm) = [];

[names,iN] = unique(names);
masses = masses(iN);
numel = numel(iN,:);

[masses,isort] = sort(masses);
names = names(isort);
numel = numel(isort,:);

length(names)

list.names = names;
list.mass = masses;
list.numel = numel;

save("potentialListAMS.mat","list")
%%
M = 100;
iM = round(list.mass) == M;
list.names(iM)
list.mass(iM)
%%
function [seeds] = addCharge(ions,seeds)
nion = size(ions,1);
nseed = size(seeds,1);
ionmat = reshape(repmat(ions',nseed,1),[],nseed*nion)';
seeds = repmat(seeds,nion,1)+ionmat;
end