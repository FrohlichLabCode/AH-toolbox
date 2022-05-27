function dataStruct = shuffleData_FT(dataStruct)
% Function to shuffle data across trials for FieldTrip functional
% connectivity analyses

data = dataStruct.trial;
nTrials = length(data);
nChannels = size(data{1},1);
tempData = cellfun(@(x)reshape(x,1,size(x,1),size(x,2)),data,'un',0);
tempData = cell2mat(tempData');

dataShuffle = zeros(nTrials,nChannels,size(data{1},2));
for iChannel = 1:nChannels
    tempShuffle = randperm(nTrials);    
    dataShuffle(:,iChannel,:) = tempData(tempShuffle,iChannel,:);
end
for iTrial = 1:nTrials
    dataStruct.trial{iTrial} = squeeze(dataShuffle(iTrial,:,:));
end



