recPath = 'E:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\0153_AttentionTask2Ephys_01_20171020_171020_141018\toAnalyze\';
load([recPath 'lfp\lfpMat'])
load('Z:\Ferret Data\0153\BehavioralTraining\VA_training\0153_Level2b_10_20171012.mat')
ephysDir = 'E:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\0153_AttentionTask2Ephys_01_20171020_171020_141018\toAnalyze\';
cd(ephysDir); % change directory to save data
% load lfpData
load([ephysDir 'lfp\lfpMat.mat'])
% load trigger data
load([ephysDir 'triggerData.mat'])
    
%% preprocess session behav data
fs = lfpFs;
twin = [-4 1]; % [-2-nanmedian(session_output_data.BehavData(:, 10)) 1]; %<<<-------------------------------- median reaction time
onInitInd = 1;
stimTouchInd = 2;

rawFs = 30000;
trialOnset = find(diff(triggerData(onInitInd,:))==1)./rawFs;
trialInit = find(diff(triggerData(onInitInd,:))==-1)./rawFs;
stimOnset = find(diff(triggerData(stimTouchInd,:))==1)./rawFs;
touch = find(diff(triggerData(stimTouchInd,:))==-1)./rawFs;

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

hitTrials = find(session_output_data.BehavData(:,Col_HitMiss) == 1 );

event = touch(hitTrials);
 
%% plot waveform to find good channels
waveformPlot = 0;
if waveformPlot == 1
    for ichan = 1:16
        load([recPath 'spikes\spk_' num2str(ichan+40)]);

        subplot(4,4,ichan)
        shadedErrorBar(1:size(spkWav,2),mean(spkWav),std(spkWav),'b')
        axis tight
        title(['chan ' num2str(ichan)])
    end
end
%% pick a channel 
figure();
for i = 1:16
    ichan = 64 + i; % VC channel starts at 65
    load([recPath 'spikes\spk_' num2str(ichan)]);
    numBins = 50;
    binWidth = (twin(2)-twin(1))/numBins; % in seconds
    numSpk = length(spkTime);

    allSpks = [];
    for iEvt = 1:length(event)
        OnsetTime = event(iEvt);
        startTime = OnsetTime + twin(1);
        endTime = OnsetTime + twin(2);

        tmpSpks = spkTime(spkTime>startTime & spkTime<endTime);
        spkZerod = tmpSpks - OnsetTime;

        allSpks = [allSpks spkZerod];

        RasterSpikeTimes{iEvt} = spkZerod;
    end
    subplot(4,4,i)
    N = histcounts(allSpks,numBins)./binWidth./length(event);
    plot(N) %bar(N)
    hold on

    xlabel('Time [Sec]')
    ylabel('Firing Rate [Hz]')
    title(['chan ' num2str(i)]) % 'Firing Rate; Bin width: ' num2str(binWidth) 's'])
    axis tight
    set(gca,'XTick',linspace(1,numBins,6))
    set(gca,'XTickLabel',linspace(twin(1),twin(2),6))

end

