mlim = [99.5 600.5];
el = ["C" "H" "O" "N" "F"];
charge = "-";
OtoC = [0.1 2.8];
HtoC = [0.6 2];
ul = [20 36 24 3 20];
ll = [4 0 2 0 0];

masses = [];
numel = [];
names = strings;

% Group1: Fluorinated carboxylic acids and dicarboxylic acids.
% CH2O2, CHO2F, CO2F2
seeds = [1 2 2 0 0; 1 1 2 0 1; 1 0 2 0 2];
% Charged by cluster with nitrate or deprotonation (defluorization?)
ions = [0 0 3 1 0; 0 -1 0 0 0];
% Only add CF2, or at most one more acid group, at most one O or CHF
clusters = [1 0 0 0 2; 1 0 2 0 0; 1 0 0 0 0;  0 0 1 0 0; 1 1 0 0 1];
nmax = [20 1 0 1 1];

seeds = addCharge(ions,seeds);

[mass,nums,name] = buildComp(el,seeds,clusters,nmax,mlim,charge);
masses = [masses; mass];
numel = [numel; nums];
names = [names; name];

% Group2: Carboxylic and dicarboxylic acids, with double bonds or ring
% structures, and nitrate groups.
% CH4, CH2O2, C2H204
seeds = [1 4 0 0 0; 1 2 2 0 0; 2 2 4 0 0];
% Charged by cluster with nitrate or deprotonation
ions = [0 0 3 1 0; 0 -1 0 0 0];
% Addition of CH2, CO, C, O, CH2O, NO2(-H), NH, or HNO3
clusters = [1 2 0 0 0; 1 0 1 0 0; 1 0 0 0 0; 0 0 1 0 0;...
    1 2 1 0 0; 0 -1 2 1 0; 0 1 0 1 0; 0 1 3 1 0];
nmax = [20 10 3 4 20 2 2 1];

seeds = addCharge(ions,seeds);

[mass,nums,name] = buildComp(el,seeds,clusters,nmax,mlim,charge);
masses = [masses; round(mass,4)];
numel = [numel; nums];
names = [names; name];

% Group3: Carboxylic and dicarboxylic acids, with double bonds or ring
% structure. Radicals.
% CH4, CH2O2, C2H204
seeds = [1 3 0 0 0; 1 1 2 0 0; 2 1 4 0 0];
% Charged by cluster with nitrate or deprotonation
ions = [0 0 3 1 0];
% Addition of CH2, CO, C, CH2O, O
clusters = [1 2 0 0 0; 1 0 1 0 0; 1 0 0 0 0; 1 2 1 0 0; 0 0 1 0 0];
nmax = [20 5 3 20 4];

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
%%
% For testing purposes:
% iterm = mlim(1);
% namList = {};
% massList = [];
% while true
%     iterm = iterm+0.02;
%     if iterm >= mlim(2)
%         break
%     end
%     [namL,mL] = find_composition(iterm);
%     namList = [namList namL];
%     massList = [massList; round(mL,4)];
% end
% 
% namList = modNams(string(namList'));
% %names = modNams(names);
% [~,imiss] = setdiff(namList,names);
% [~,iextra] = setdiff(names,namList);
% missed = namList(imiss);
% extra = names(iextra);
% fprintf("\nFormulas not in new method: %i",length(imiss))
% fprintf("\nAdditional formulas in new method: %i",length(iextra))
% fprintf("\nTotal number of formulas in new method: %i",length(names))
% fprintf("\nTotal number of formulas in old method: %i",length(namList))
% fprintf("\n")
%%
potList = struct;
potList.names = names;
potList.mass = masses;
potList.numel = numel;
save("potentialList100_600.mat","potList")
%%
% fitList = modNams(list.compNam);
% %fitList = string(list.compNam);
% iuse = list.peakList > mlim(1) & list.peakList < mlim(2);
% %size(fitList(iuse))
% n1 = setdiff(fitList(iuse),names)
% n2 = setdiff(fitList(iuse),namList)
% n3 = setdiff(n1,n2)
%%
function [seeds] = addCharge(ions,seeds)
nion = size(ions,1);
nseed = size(seeds,1);
ionmat = reshape(repmat(ions',nseed,1),[],nseed*nion)';
seeds = repmat(seeds,nion,1)+ionmat;
end