function [meanPower,maxFreq] = AH_timeWindowMean(freqBands,timeWins,foi,tvec,powerMat)
% This function takes in a matrix of value (eg.power or PLV) and calculate
% mean across freqBand and timeWin
% Input: eg. freqBands = {[4,8],[10,14],[32,70]}; % theta, alpha, gamma
%            timeWins  = {[-5,-4],[-3,0],[1,2]}; % before, during, after opto
%
% Created by AH 2/2022

meanPower = nan(numel(freqBands),numel(timeWins));
maxFreq = nan(numel(freqBands),numel(timeWins));

% if funcCon doesn't have foi saved, use this, but double check if it
% matches original frequencies. If funcCon has foi, use funcCon.foi
% numFreqs = 100;
% lowFreq  = 2;
% highFreq = 30;
% foi      = linspace(lowFreq,highFreq,numFreqs); % linear spacing

% numFreqs = 150;
% lowFreq  = 2;
% highFreq = 128;
% foi      = logspace(lowFreq,highFreq,numFreqs); % log spacing


for iFreq = 1:size(freqBands,1)
    for iTwin = 1:size(timeWins,1)
        fMask = foi>freqBands{iFreq}(1) & foi<=freqBands{iFreq}(2);
        tMask = tvec>timeWins{iTwin}(1) & tvec<=timeWins{iTwin}(2);
        xslice = powerMat(fMask,tMask);
        meanPower(iFreq,iTwin) = mean(mean(xslice,2),1); %average over time, sum over freq
        
        % Get the freq with the max power
        [~,maxID] = max(xslice(:)); %[max,maxID]
        [X,~]=ind2sub(size(xslice),maxID); %[X,Y]
        foi_band = foi(fMask); % get freq of interest within freq band
        maxFreq(iFreq,iTwin) = foi_band(X); % get max freq
    end
end
end


