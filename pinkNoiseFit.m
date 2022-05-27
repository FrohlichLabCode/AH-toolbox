function [powerSlope,powerYintercept,pinkNoise] = pinkNoiseFit(frequencyDomain,channelByPower,varargin)
%% Estimate the power slope for the input data
% INPUTS
%   frequencyDomain - [1 x numFrequencies] vector of frequencies
%   channelByPower - [numChannels x numFreq], output from Fouier transform
%   options - structure of options
%
%   OPTIONS
%       cutAlphaPower  - 0/1, better estimate without the alpha band included
%           default to cut the alpha power out of the estimate
%       downSample - maybe get rid of this
%
% OUTPUTS
%   pinkSlopeFit - [1 x numChannels] slope by channel
%   pinkYintercept - [1 x numChannels] y-intercept from fit
%   pinkNoise - [numChannels x numFreq] estimate of the background noise by frequency
%
% To get denoised spectrum:
% newSpec = spec - pinkNoise;


% Get number of channels and number of frequencies from channelByPower input
[numChannels,numFrequencies] = size(channelByPower);

if ~isempty(varargin)
    options = varargin{1};
else
    options = struct([]);
end
if isfield(options,'cutAlphaPower')
    cutAlphaPower = options.cutAlphaPower;
    % choose 0 if you do not want to cut the alpha power out
else
    % default is the canonical 8 to 12 hertz range
    cutAlphaPower = [10 18]; % for ferret
    % recommendation: cut from 4 to 30 because beta and theta oscillations are common
end

% Down sample the power
% go to half frequency increments
if isfield(options,'downSampleFreq')
    downSampleFreq = options.downSampleFreq;
else
    downSampleFreq = 0;
end

% prepare output variables
powerSlope = NaN(1,numChannels);
powerYintercept = NaN(1,numChannels);
pinkNoise = NaN(numChannels,numFrequencies);

% fit a line to the log log - this is the 1/f slope for the subject
for chanIdx = 1:numChannels
    
    %% Plot frequency by log power for each channel
    power = channelByPower(chanIdx,:);
    
     %% Down sample the data - using mean of every 1/2 frequency window
    if downSampleFreq
        numFreqPoints = 1000;
        [est_frequencyDomain,est_power] = downSampleFreq(frequencyDomain,power,numFreqPoints,frequencyDomain(1),frequencyDomain(end));
    else
        est_frequencyDomain = frequencyDomain;
        est_power = power;
    end
    
    %% Band cut the alpha frequency out of our estimate
    if cutAlphaPower
        alphaLow = cutAlphaPower(1);
        alphaHigh = cutAlphaPower(2);
        alphaLowIdx = find(est_frequencyDomain>=alphaLow,1);
        alphaHighIdx = find(est_frequencyDomain>alphaHigh,1);
        cutFreq = est_frequencyDomain;
        cutFreq(alphaLowIdx:alphaHighIdx) = [];
        cutPower = est_power;
        cutPower(alphaLowIdx:alphaHighIdx)= [];
    else
        % Advise to run the cut...
        % The bump in alpha power is so strong that it throws off the estimate
        cutFreq = est_frequencyDomain;
        cutPower = est_power;
    end
    
    %% run a linear fit on the cut frequency by log power
    specFit = polyfit(log(cutFreq),log(cutPower),1);
    
    % output the y-intercept and slope
    powerSlope(chanIdx) = specFit(1);
    powerYintercept(chanIdx) = specFit(2);
    
    % generate the log power values
    channelPinkNoise = exp(polyval(specFit,log(frequencyDomain)));
    
    % place back in standard values
    pinkNoise(chanIdx,:) = channelPinkNoise;
end
end