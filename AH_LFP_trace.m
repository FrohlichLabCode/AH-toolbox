% Get LFP traces for 0153
load('E:\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\Ephys\0153\toAnalyze\0153_AttentionTask3_13\lfp\lfpMat.mat')
load('E:\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\Ephys\0153\toAnalyze\0153_AttentionTask3_13\triggerData.mat')

fs = lfpFs;
onInitInd = 1;
stimTouchInd = 2;

rawFs = 30000;
trialOnset = find(diff(triggerData(onInitInd,:))==1)./rawFs;
trialInit = find(diff(triggerData(onInitInd,:))==-1)./rawFs;
stimOnset = find(diff(triggerData(stimTouchInd,:))==1)./rawFs;
touch = find(diff(triggerData(stimTouchInd,:))==-1)./rawFs;

region1Chn = 1:32;%9:24;  %<<<------------ PPC
region2Chn = 33:48; %40:56; %33:64; 33:48  %<<<------------ LP/Pulvinar
region3Chn = 49:64; % VC
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
event = event(2:end);
twin = [-4 2];
tsamps = round(twin*fs);
iev = 2; % trial number
evSamp = round(event(iev)*fs);
tvec = [twin(1):1/fs:twin(2)];
% this may not work if the window looks for samples out of range
evLFP = lfpMat(:,evSamp+tsamps(1):evSamp+tsamps(2));
figure()
subplot(3,1,1)
plot(tvec, evLFP(1,:));ylabel('Amplitude [uV]');
title('PPC');ylim([-25,25]);

subplot(3,1,2)
plot(tvec, evLFP(50,:));ylabel('Amplitude [uV]');
title('V1');ylim([-80,70]);

subplot(3,1,3)
plot(tvec, evLFP(36,:));ylabel('Amplitude [uV]');
title('LPl');ylim([-140,160]);
xlabel('Time to touch [sec]');
