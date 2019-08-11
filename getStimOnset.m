function onset = getStimOnset(triggerData, photoRaw, photoBaseline)
% This script is used for visual presentation stimulus that has triggerData as digital
% input and optional photodiode data as analytical inputto INTAN system

% measured from 2P room computer: on average, trigger data is 0.0529 sec earlier than visual stimuli onset
defaultOffset = 0.05; 

rawFs = 30000; % default sampling frequency from INTAN system
ttlOnset  = find(diff(triggerData==1)./rawFs; % in sec
ttlOffset = find(diff(triggerData==-1)./rawFs;

if nargin == 2
    photoBaseline = 0.12; % default threshold when using 0V ground; photoThreshold = 3.2; % used 3.3V ground
end 

if nargin == 1 % only has triggerData, then use triggerData to approximate
    onset = ttlOnset + defaultOffset; 
else
    photoOnset = [];
    for i = ttlOnset*rawFs % for each trigger time point
        rawAvg = movmean(photoRaw(i:i+rawFs*0.1), 200);% about every <200 sample, photodiode has a dip, this will smooth it out
        offset = find(rawAvg < photoBaseline + 0.04,1,'first'); %moving average, use 0.16 as threshold %0169 use:find(diff(temp(i:i+rawFs*0.1))==1,1,'last');
        if isempty(offset); offset = 0;end
        iphotoOnset = (i+offset)/rawFs; 
        photoOnset = [photoOnset, iphotoOnset]; % the last drop->0 is when gray->black in photodiode
    end
    onset = photoOnset;
end

