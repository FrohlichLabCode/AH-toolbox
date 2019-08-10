function [regionPairNames, regionPair_Names, regionPairNamesGC] = getRegionPairName(regionNames, regionPairs)

% generate name pair cells
for iPair = 1:numel(regionPairs)
    regionPairNames{iPair} = [regionNames{regionPairs{iPair}(1)} '-' regionNames{regionPairs{iPair}(2)}];
    regionPair_Names{iPair} = [regionNames{regionPairs{iPair}(1)} '_' regionNames{regionPairs{iPair}(2)}];%{'FC_PPC', 'LPl_PPC', 'LPl_VC', 'PPC_VC'};
    regionPairNamesGC{2*iPair-1} = [regionNames{regionPairs{iPair}(1)} '->' regionNames{regionPairs{iPair}(2)}];%{'LPl->PPC','PPC->LPl';'LPl->VC','VC->LPl'};
    regionPairNamesGC{2*iPair} = [regionNames{regionPairs{iPair}(2)} '->' regionNames{regionPairs{iPair}(1)}];
end
return
end