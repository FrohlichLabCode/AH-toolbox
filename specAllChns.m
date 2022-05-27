function  [spec, specgram, f, t] = specAllChns(lfpMat, params)
% This funcion takes nChn x nT matrix of lfpMat and calculate spectrogram
% for each channel, using a sliding window, then average time to get spec
% AH 2020/9
if nargin == 2
    window   = params.window; 
    noverlap = params.noverlap;
    foi      = params.foi; 
    lfpFs    = params.lfpFs;
    numFreqs = params.numFreqs;
else % default values
    window   = 5*1024; % ~5s 
    noverlap = 0.5*window;
    numFreqs = 100;
    lowFreq  = 2;
    highFreq = 80;
    %foi      = logspace(log10(lowFreq),log10(highFreq),numFreqs);
    foi      = linspace(lowFreq, highFreq, numFreqs);
    lfpFs    = 1000;
end
% Initiate empty arrays
spec     = nan(size(lfpMat,1),numFreqs);
[s,f,t]  = spectrogram(lfpMat(1,:),window,noverlap,foi,lfpFs); % just to get t, f values
specgram = nan(size(lfpMat,1), size(s,1), size(s,2));

parfor ichan = 1:size(lfpMat,1)
    %display(['computing spec ' num2str(ichan)])
    [s1,~,~] = spectrogram(lfpMat(ichan,:),window,noverlap,foi,lfpFs);
    spec(ichan,:) = mean(abs(s1).^2,2); % LFP power
    specgram(ichan, :, :) = s1;
end
end