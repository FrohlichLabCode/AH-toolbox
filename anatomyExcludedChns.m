function [exLFPChn, exSUChn] = anatomyExcludedChns(animalCode)
% number coming from excel sheet ChannelMap 2020/10/26
for iRegion = 1:4
    exLFPChn{iRegion} = [];
    exSUChn{iRegion} = [];
end
switch animalCode
    case '0171'
        exLFPChn{2} = [1,2];        
        exSUChn{2} = [5,6];
    case '0179'
        exLFPChn{2} = [1,2];        
        exSUChn{2} = [5,6];
    case '0180'
        exLFPChn{2} = [1,2,8,16];
        exSUChn{2} = [5,6,12,13];
    case '0181' % 10 channels excluded;
        exLFPChn{2} = [1,6:11,14:16];
        exSUChn{2} = [2:5,10:15];        
end