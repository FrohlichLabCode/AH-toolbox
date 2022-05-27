%% Get LFP traces for 0153
% load('C:\Users\angel\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\0153_AttentionTask2Ephys_01_20171020_171020_141018\toAnalyze\lfp\lfpMat.mat')
% load('C:\Users\angel\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\0153_AttentionTask2Ephys_01_20171020_171020_141018\toAnalyze\triggerData.mat')
load('E:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\toAnalyze\0153_AttentionTask2Ephys_01_20171020_171020_141018\lfp\lfpMat.mat')
load('E:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\toAnalyze\0153_AttentionTask2Ephys_01_20171020_171020_141018\triggerData.mat')

fs = lfpFs;
onInitInd = 1;
stimTouchInd = 2;

rawFs = 30000;
trialOnset = find(diff(triggerData(onInitInd,:))==1)./rawFs;
trialInit = find(diff(triggerData(onInitInd,:))==-1)./rawFs;
stimOnset = find(diff(triggerData(stimTouchInd,:))==1)./rawFs;
touch = find(diff(triggerData(stimTouchInd,:))==-1)./rawFs;

region1Chn = 1:32;%9:24;  %<<<------------ PPC
region2Chn = 41:56; %40:56; %33:64; 33:48  %<<<------------ LP/Pulvinar
region3Chn = 65:80; % VC
numChnReg1 = numel(region1Chn);
numChnReg2 = numel(region2Chn);
numChnReg3 = numel(region3Chn);

region1LFP = lfpMat(region1Chn,:);
region2LFP = lfpMat(region2Chn,:);
region3LFP = lfpMat(region3Chn,:);

% behav data column names
Col_TrialNumber = 1;
Col_StimulusWindow = 2;
Col_TrialOnsetToInit = 4; % How long the spout light is on prior to trial initiation
Col_XCoordinateTouch = 5;
Col_YCoordinateTouch = 6;
Col_TouchTimeStampHr = 7;
Col_TouchTimeStampMin = 8;
Col_TouchTimeStampSec = 9;
Col_TrialOnsetToTouch = 10;
Col_HitMiss = 11;   % 0 for miss or 1 for touch; 2 for premature touch

event = touch;
twin = [-5 0];
tsamps = round(twin*fs);

% define a band pass filter at 4-6Hz
%downFac = 30;         % factor to downsample raw signal for LFP
bpFreq = [4,6];       % lowpass frequency
ford    = 4;          % filter order
%lfpFs   = Fs/downFac; % LFP sample rate
nyqFreq = fs/2;       % Nyquest frequency of raw signal
[a,b]   = butter(ford,bpFreq/nyqFreq,'bandpass'); % lowpass filter 
%[dr,nr] =  rat(Fs/(Fs/downFac)); % Ratio for downsampling
peakError = cell(length(region2Chn),45);

for ichn = region2Chn
    for iev = 1:45  % trial number
    tic
    evSamp = round(event(iev)*fs);
    tvec = [twin(1):1/fs:twin(2)];
    % this may not work if the window looks for samples out of range
    evLFP = lfpMat(:,evSamp+tsamps(1):evSamp+tsamps(2));

    %% plot event window LFP for 3 areas
    % figure()
    % subplot(3,1,1)
    % plot(tvec, evLFP(1,:));ylabel('Amplitude [uV]');
    % title('PPC');ylim([-300,300]);
    % subplot(3,1,2)
    % plot(tvec, evLFP(41,:));ylabel('Amplitude [uV]');
    % title('Pulvinar');ylim([-300,300]);
    % subplot(3,1,3)
    % plot(tvec, evLFP(65,:));ylabel('Amplitude [uV]');xlabel('Time to touch [sec]');
    % title('VC');ylim([-300,300]);

    %% Peak detection
    
    % filter LFP signal
    evLFPTheta = filtfilt(a,b,evLFP(ichn,:));
%     figure()
%     plot(tvec, evLFPTheta); % plot 1 channel filtered data -- looks smooth

    % detect peak
    temp = diff(evLFPTheta);
    temp(temp>0) = 1; % value increasing
    temp(temp<0) = -1; % value decreasing
    % figure() % temp from 1 to -1 is the peak location
    % plot(tvec(1:length(temp)),temp);
    myPeakInds = find((diff(temp)==-2));
    myPeaks = zeros(size(evLFPTheta));
    myPeaks(myPeakInds)=abs(evLFPTheta(myPeakInds));
%     figure() % temp from 1 to -1 is the peak location
%     plot(tvec,myPeaks);    

    % % calculate peak based on every 3 peaks moving window (600 ms)  --
    % % somehow skip 1 spike
    % window = 0.6*fs; % 0.6s time window
    % nWin = 5*(twin(2)-twin(1))-3; 
    % nextPeakInds = [];
    % for iWin = 1:nWin
    %     t1 = fs*(iWin-1)/5+1;
    %     t2 = fs*(iWin-1)/5+window;
    %     prevPeakInd = find(myPeaks(t1:t2));
    %     p2p = round(mean(diff(prevPeakInd))); % peak-to-peak interval
    %     nextPeakInd = prevPeakInd(length(prevPeakInd)) + p2p; % take the last peak time + p2p
    %     nextPeakInds(iWin) = nextPeakInd+t1;
    % end

    % calculate peak based on every 3 peaks moving window (600 ms)
    window = 3; % 0.6s time window
    nWin = length(myPeakInds) - window + 1; 
    nextPeakInds = [];
    for iWin = 1:nWin
        t1 = iWin;
        t2 = iWin + window -1;
        p2p = round(mean(diff(myPeakInds(t1:t2)))); % peak-to-peak interval
        nextPeakInds(iWin) = myPeakInds(t2) + p2p; % take the last peak time + p2p
    end

    nextPeakInds = nextPeakInds(1:(length(nextPeakInds)-1)); % exclude the last peak (out of range)
    nextPeaks = nan(size(evLFPTheta));
    nextPeaks(nextPeakInds)= 1;
    peakError{ichn-40,iev} = myPeakInds(4:length(myPeakInds))-nextPeakInds;
    tocTrial = toc;
%     figure(); % overlap original signal, real peak with predicted peak
%     hold on
%     plot(tvec, evLFP(41,:));
%     myPeaks(myPeaks>0) = 1;
%     myPeaks(myPeaks==0) = NaN;
%     scatter(tvec, myPeaks*175, 'filled','d');
%     scatter(tvec, nextPeaks*150, 'filled'); % move dots to 150uV
%     legend('LFP','real peak','predicted peak');
%     xlabel('Time to touch [sec]');
%     ylabel('LFP [uV]');
%     ylim([-200,200]);
%     box on;
    end
end
    % calculate difference between 2 array

% convert all cell into mat
peakErrorMat = [];
for ichn = 1:length(region2Chn)
    for iev = 1:45  % trial number
        peakErrorMat = [peakErrorMat,peakError{ichn,iev}];
    end
end
figure();
histogram(abs(peakErrorMat)/200 * 100,'BinWidth',2,'normalization','probability');
xlabel('Peak time difference (%)');
ylabel('% samples');
xlim([0,40]);

% 24.09+24.85+17.81+10.96+7.93
