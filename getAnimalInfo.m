function region = getAnimalInfo(animalCode)
% This function takes in animalCode in string and outputs regionInfo struct
% Created by AH 2021

switch animalCode
    case {'0147'}
        region.Names   = {'LPl','PPC','VC','EEG'};
        region.IDs     = [1,2,3,4];
        region.N       = numel(region.Names);
        region.allChns = {[65:80],[1:32],[33:64],[7:8] + 80};
        region.PairIDs = {[1,2],[1,3],[2,3],[1,4],[2,4],[3,4]};
        region.NPair   = numel(region.PairIDs);
        [region.PairNames, region.Pair_Names, region.PairNamesGC, region.Pair_NamesGC] = getRegionPairName(region.Names, region.PairIDs);
    case {'0139','0151','0153'}
        region.Names   = {'LPl','PPC','VC','EEG'};
        region.IDs     = [1,2,3,4];
        region.N       = numel(region.Names);
        region.allChns = {[49:64],[1:32],[33:48],[7:8] + 64};
        region.PairIDs = {[1,2],[1,3],[2,3],[1,4],[2,4],[3,4]};
        region.NPair   = numel(region.PairIDs);
        [region.PairNames, region.Pair_Names, region.PairNamesGC, region.Pair_NamesGC] = getRegionPairName(region.Names, region.PairIDs);
    case {'0158'}
        region.Names   = {'LPl','PPC','VC','EEG'};
        region.IDs     = [1,2,3,4];
        region.N       = numel(region.Names);
        region.allChns = {[33:48],...
            [17 19 20 22 23 25 27 29 30 31 32],...
            [1 2 4 5 7 8 9 10 11 15 16],...
            [7:8] + 48};
        %region.allChns = {[33:48],[16:32],[1:16],[7:8] + 48};
        region.PairIDs = {[1,2],[1,3],[2,3],[1,4],[2,4],[3,4]};
        region.NPair   = numel(region.PairIDs);
        [region.PairNames, region.Pair_Names, region.PairNamesGC, region.Pair_NamesGC] = getRegionPairName(region.Names, region.PairIDs);

    case {'0181','0180','0179','0171','0172','0173'}
        region.Names   = {'PFC','LPl','PPC','VC'};
        region.IDs     = [1,2,3,4];
        region.N       = numel(region.Names);
        region.allChns = {[1:16],[17:32],[33:48],[49:64]};
        region.PairIDs = {[2,3],[1,2],[1,3],[1,4],[2,4],[3,4]};
        region.NPair   = numel(region.PairIDs);
        [region.PairNames, region.Pair_Names, region.PairNamesGC, region.Pair_NamesGC] = getRegionPairName(region.Names, region.PairIDs);
    case {'0182','0183','0185'} % alpha opto in PPc, arnoldTongue control
        region.Names   = {'PPC'};
        region.IDs     = [1];
        region.N       = numel(region.Names);
        region.allChns = {[1:16]};
%         region.PairIDs = {[]};
%         region.NPair   = numel(region.PairIDs);
%         [region.PairNames, region.Pair_Names, region.PairNamesGC, region.Pair_NamesGC] = getRegionPairName(region.Names, region.PairIDs);
    case {'0186','0187','0188'} % arnoldTongue control, bonescrew tACS, tACS-opto
        region.Names   = {'PPC','VC'};
        region.IDs     = [1,2];
        region.N       = numel(region.Names);
        region.allChns = {[1:16],[17:32]};
        region.PairIDs = {[1,2]};
        region.NPair   = numel(region.PairIDs);
        [region.PairNames, region.Pair_Names, region.PairNamesGC, region.Pair_NamesGC] = getRegionPairName(region.Names, region.PairIDs);
    case {'0201','0202','0200','0203'} 
        region.Names   = {'PMC','CLA','PPC'};
        region.IDs     = [1,2,3];
        region.N       = numel(region.Names);
        region.allChns = {[1:16],[17:32],[33:48]};
        region.PairIDs = {[1,2],[1,3],[2,3]};
        region.NPair   = numel(region.PairIDs);
%         switch animalCode
%                 case {'0201'}  
%                     %if (strcmp(level,'7') && ismember(sessionID,{'09','11','13'}))
%                     if (strcmp(level,'7') && ismember(sessionID,{'09'}))
%                         region.Names = {'PMC','CLA'};
%                         region.IDs     = [1,2];
%                         region.N       = numel(region.Names);
%                         region.allChns = {[1:16],[17:32]};
%                         region.PairIDs = {[1,2]};
%                         region.NPair   = numel(region.PairIDs);
% %                     end
%                 case {'0200'} 
%                     if (strcmp(level,'7') && ismember(sessionID,{'01','06','07','08'}))
%                         region.Names = {'PMC','CLA'};
%                         region.IDs     = [1,2];
%                         region.N       = numel(region.Names);
%                         region.allChns = {[1:16],[17:32]};
%                         region.PairIDs = {[1,2]};
%                         region.NPair   = numel(region.PairIDs);
%                     end
%                 case {'0203'} 
%                     %if (strcmp(level,'7') && ismember(sessionID,{'02','03','05','06','07','08','09','11'}))
%                     if (strcmp(level,'7') && ismember(sessionID,{'07'}))
%                         region.Names = {'PMC','CLA'};
%                         region.IDs     = [1,2];
%                         region.N       = numel(region.Names);
%                         region.allChns = {[1:16],[17:32]};
%                         region.PairIDs = {[1,2]};
%                       region.NPair   = numel(region.PairIDs);
%                   end
%         end
            end
        [region.PairNames, region.Pair_Names, region.PairNamesGC, region.Pair_NamesGC] = getRegionPairName(region.Names, region.PairIDs);

end
