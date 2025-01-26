nlists = 5;
path = "/home/vami/Documents/MATLAB/PeakFitAlgorithm_2023/potentialCompounds/testingLists/";
prefix = "testList";
mlim = [199.5 400.5];
el = ["C" "H" "O" "N" "F" "X"];
charge = "-";
OtoC = [0.1 2.8];
HtoC = [0.6 2];
ul = [20 30 20 3 20 nlists];
ll = [4 0 2 0 0 0];
for i = 1:nlists
    masses = [];
    numel = [];
    names = strings;

    % Group1: Fluorinated carboxylic acids and dicarboxylic acids.
    % CH2O2, CHO2F, CO2F2
    seeds = [1 2 2 0 0 0; 1 1 2 0 1 0; 1 0 2 0 2 0];
    % Charged by cluster with nitrate or deprotonation (defluorization?)
    ions = [0 0 3 1 0 0; 0 -1 0 0 0 0; 0 0 0 0 -1 0];
    % Only add CF2, or at most one more acid group, at most one O or CHF
    clusters = [1 0 0 0 2 0; 1 0 2 0 0 0; 1 0 0 0 0 0;  0 0 1 0 0 0; 1 1 0 0 1 0; 0 0 0 0 0 1];
    nmax = [20 1 0 1 1 0];

    nmax(end) = i-1;

    seeds = addCharge(ions,seeds);

    [mass,nums,name] = buildComp(el,seeds,clusters,nmax,mlim,charge);
    masses = [masses; mass];
    numel = [numel; nums];
    names = [names; name];

    % Group2: Carboxylic and dicarboxylic acids, with double bonds or ring
    % structures, and nitrate groups.
    % CH4, CH2O2,
    seeds = [1 2 2 0 0 0; 2 2 4 0 0 0];
    % Charged by cluster with nitrate or deprotonation
    ions = [0 0 3 1 0 0; 0 -1 0 0 0 0];
    % Addition of CH2, CO, C, O, NO2(-H), NH, HNO3 or H
    clusters = [1 2 0 0 0 0; 1 0 1 0 0 0; 1 0 0 0 0 0; 0 0 1 0 0 0;...
        0 -1 2 1 0 0; 0 1 0 1 0 0; 0 1 3 1 0 0; 0 1 0 0 0 0; 0 0 0 0 0 1];
    nmax = [20 5 3 12 2 2 1 0 0];

    nmax(end) = i-1;

    seeds = addCharge(ions,seeds);

    [mass,nums,name] = buildComp(el,seeds,clusters,nmax,mlim,charge);
    masses = [masses; round(mass,4)];
    numel = [numel; nums];
    names = [names; name];

    % Group3: Carboxylic and dicarboxylic acids, with double bonds or ring
    % structure. Radicals.
    % CH4, CH2O2.
    seeds = [1 1 2 0 0 0; 2 1 4 0 0 0];
    % Charged by cluster with nitrate or deprotonation
    ions = [0 0 3 1 0 0];
    % Addition of CH2, CO, C, O
    clusters = [1 2 0 0 0 0; 1 0 1 0 0 0; 1 0 0 0 0 0; 0 0 1 0 0 0; 0 1 0 0 0 0; 0 0 0 0 0 1];
    nmax = [20 5 3 12 0 0];

    nmax(end) = i-1;

    seeds = addCharge(ions,seeds);

    [mass,nums,name] = buildComp(el,seeds,clusters,nmax,mlim,charge);
    masses = [masses; round(mass,4)];
    numel = [numel; nums];
    names = [names; name];

    names = names(2:end);
    OCratio = (numel(:,3)-(max(numel(:,4),1)-1)*3)./numel(:,1);
    HCratio = (numel(:,2)+numel(:,5))./numel(:,1);
    irm = OCratio < OtoC(1) | OCratio > OtoC(2);
    irm = irm | HCratio < HtoC(1) | HCratio > HtoC(2);
    irm = irm | any(numel>ul,2);
    irm = irm | any(numel<ll,2);
    irm = irm | (numel(:,4)==3 & rem(numel(:,2),2)==1);
    irm = irm | (numel(:,5)==0 & numel(:,3)<4);

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

    potList = struct;
    potList.names = names;
    potList.mass = masses;
    potList.numel = numel;
    fnam = sprintf("%s%s%i.mat",path,prefix,i-1);
    save(fnam,"potList")

end
%%
function [seeds] = addCharge(ions,seeds)
nion = size(ions,1);
nseed = size(seeds,1);
ionmat = reshape(repmat(ions',nseed,1),[],nseed*nion)';
seeds = repmat(seeds,nion,1)+ionmat;
end