function region = getAnimalInfo(animalCode)
% This function takes in animalCode in string and outputs regionInfo struct
switch animalCode
    case {'0171','0172','0173'}
        region.Names   = {'PFC','LPl','PPC','VC'};
        region.IDs     = [1,2,3,4];
        region.N       = numel(region.Names);
        region.allChns = {[1:16],[17:32],[33:48],[49:64]};
        region.PairIDs = {[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]};
        region.NPair   = numel(region.PairIDs);
        [region.PairNames, region.Pair_Names, region.PairNamesGC] = getRegionPairName(region.Names, region.PairIDs);
end
