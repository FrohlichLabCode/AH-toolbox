clear
tic

cluster = 0;
skipRec = 1;
linORlog = 2; %freqs of interest: 1=linear 2=log
MedianorPCA = 3; 
animalCode = '0147';


if cluster == 0 %linux use '/', windows matlab can use both '/' and '\'
    addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys'));
    PreprocessDir = ['E:/FerretData/' animalCode '/Preprocessed/'];
    AnalysisDir   = ['E:/FerretData/' animalCode '/Analyzed/'];
    BehavDatDir   = ['E:/FerretData/' animalCode '/behav/'];
    GroupAnalysisDir = ['E:/FerretData/' animalCode '/GroupAnalysis/Power/'];
elseif cluster == 1
    addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
    PreprocessDir = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Preprocessed/']; % CHANGE FOR KILLDEVIL VS LONGLEAF
    AnalysisDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Analyzed/'];
    BehavDatDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/behav/'];

%     %code for initialising parallel computing
%     numCore = 8; % USR DEFINE
%     parpool('local',numCore);  
end

% Different folder name if use median or PCA LFP
if MedianorPCA == 0;     folderSuffix = '_validChns_150f'; %use all valid channels
% elseif MedianorPCA == 1; folderSuffix = '_median'; %use median of all valid channels
% elseif MedianorPCA == 2; folderSuffix = '_PCA';
elseif MedianorPCA == 3; folderSuffix = '_random'; % randomly pick 1 channel
end

fileInfo = dir([PreprocessDir animalCode '_AttentionTask6*']); % detect files to load/convert  '_LateralVideo*'

% loop through each recording
for irec = 1:numel(fileInfo)
    recName = fileInfo(irec).name;
    splitName   = strsplit(recName,'_');
    sessionName = [splitName{2}(14:end) splitName{3}];
    %if cluster == 0 && datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180726', 'InputFormat', 'yyyyMMdd'); continue;end

    rootPreprocessDir = [PreprocessDir recName '/'];
    rootAnalysisDir   = [AnalysisDir recName '/FC' folderSuffix '/']; %eg. FC_median
    % below 3 lines are dealing with multiple record with the same name,
    % use the first record
    BehavNames = dir([BehavDatDir splitName{1} '_Level' splitName{2}(14:end) '_' splitName{3} '_' splitName{4} '*']);
    BehavName  = BehavNames.name;
    rootBehavDatDir   = [BehavDatDir BehavName];

%     if ~exist(join(rootAnalysisDir),'dir') 
%         mkdir(join(rootAnalysisDir)); fprintf('\nWorking on record %s =============== \n',recName'); end

    % load and process ttl data in ephys
    lfpMat = is_load([rootPreprocessDir 'lfp/lfpMat'],'lfpMat'); % don't feed in denoised data with NaN values
    lfpFs  = is_load([rootPreprocessDir 'lfp/lfpMat'],'lfpFs');
    load([rootPreprocessDir 'triggerData']);
    load(rootBehavDatDir)


%%
% Define frequencies of interest. Linear spacing for Phase slope index, and spacing for all other methods.
if linORlog == 1
    numFreqs = 100;
    lowFreq  = 1;
    highFreq = 80;
    %foi      = linspace(lowFreq,highFreq,numFreqs); % linear spacing
elseif linORlog == 2
    numFreqs = 150; %<<<--------USERDEFINE
    lowFreq  = 2;
    highFreq = 128;
    %foi   = logspace(log10(lowFreq),log10(highFreq),numFreqs); % log spacing
end

% region info
if strcmp(splitName{4}(1:8),'20170824') || strcmp(splitName{4}(1:8),'20170825')% date info
    lfp.validChn = {9:24, 33:64};
    regionNames = {'LPl','PPC'};
    regionPairs = {[1,2]};
    onInitInd = 2;
    stimTouchInd = 3;    
elseif datetime(splitName{4}(1:8), 'InputFormat', 'yyyyMMdd') >= datetime('20170829', 'InputFormat', 'yyyyMMdd') &&...
       datetime(splitName{4}(1:8), 'InputFormat', 'yyyyMMdd') <= datetime('20170901', 'InputFormat', 'yyyyMMdd');
    lfp.validChn = {1:32, 41:56};
    regionNames = {'PPC','LPl'};
    regionPairs = {[2,1]};
    onInitInd = 2;
    stimTouchInd = 3;
elseif datetime(splitName{4}(1:8), 'InputFormat', 'yyyyMMdd') >= datetime('20170908', 'InputFormat', 'yyyyMMdd') &&...
       datetime(splitName{4}(1:8), 'InputFormat', 'yyyyMMdd') <= datetime('20170922', 'InputFormat', 'yyyyMMdd');
    lfp.validChn = {1:32,33:48};
    regionNames = {'PPC','LPl'};
    regionPairs = {[2,1]};
    onInitInd = 1;
    stimTouchInd = 2;
elseif datetime(splitName{4}(1:8), 'InputFormat', 'yyyyMMdd') >= datetime('20170925', 'InputFormat', 'yyyyMMdd') &&...
       datetime(splitName{4}(1:8), 'InputFormat', 'yyyyMMdd') <= datetime('20171025', 'InputFormat', 'yyyyMMdd');
    lfp.validChn = {1:32, 33:48};
    regionNames = {'VC','LPl'};
    regionPairs = {[2,1]};
    onInitInd = 1;
    stimTouchInd = 2;
end
if MedianorPCA == 3
    lfp.validChn = {lfp.validChn{1}(1),lfp.validChn{2}(1)}; %only pick the first channel in each region
end

% get valid channels
for i = 1:numel(lfp.validChn) 
    regionChn{i} = lfp.validChn{i}; % Pulvinar, PPC, VC
    regionLFP{i} = lfpMat(lfp.validChn{i},:); % reordered channel correspond to reordered lfp
end


rawFs = 30000;
trialOnset = find(diff(triggerData(onInitInd,:))==1)./rawFs;
trialInit = find(diff(triggerData(onInitInd,:))==-1)./rawFs;
stimOnset = find(diff(triggerData(stimTouchInd,:))==1)./rawFs;
touch = find(diff(triggerData(stimTouchInd,:))==-1)./rawFs;


%% preprocess session behav data

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

hitTrials = find( session_output_data.BehavData(:,Col_HitMiss) == 1 );

condNames = {'Init','StimOnset','Touch'};
condID    = [1,3];%[1,2,3];
numConds  = numel(condID);

%% within region 
for iCond = 1:numConds
    condName = condNames{iCond};
    switch condName
        case 'Init'
            evtTime = trialInit(hitTrials);
            twin  = [-5 10]; %[-2,8]
            baseTwin = [-3 -2];
        case 'StimOnset'
            evtTime = stimOnset(hitTrials);
            twin  = [-7 5];
            baseTwin = [-7 -5.5];
        case 'Touch'
            evtTime = touch(hitTrials);
            twin  = [-4 5];
            baseTwin = [-9 -8];           
    end
    
    evtTime(evtTime < abs(twin(1)) | evtTime > size(lfpMat,2)/lfpFs-twin(2)) = []; % exclude events beyond analysis time window
    evtTimes{iCond} = evtTime;
    twins{iCond}    = twin;
    baseTwins{iCond}= baseTwin;
end

% only need to save this for the first time running
save([rootPreprocessDir 'eventTimes'], 'evtTimes','twins','baseTwins', 'condNames', 'condID', 'lfpFs');

% % for each region pairs
% regionPair_FunConn_V2(skipRec, linORlog, lowFreq, highFreq, numFreqs, lfpFs,...
%     evtTimes,twins,baseTwins, condNames, condID, regionPairs, regionNames, ...
%     regionLFP, regionChn, rootAnalysisDir)
%     %V2 allows diff twin and baseTwin

% for each region
region_spec_by_trial(skipRec, linORlog, lowFreq, highFreq, numFreqs, lfpFs,...
    evtTimes,twins,baseTwins, condNames, condID, regionNames, ...
    regionLFP, regionChn, sessionName, GroupAnalysisDir);
    
end % end of all records for an animal

x = toc;
fprintf('time required =%f sec\n', x);