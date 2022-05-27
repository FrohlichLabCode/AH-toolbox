% globalPath = 'J:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\Analyzed\';
% recs       = dir([globalPath '0153_AttentionTask3*\1-80Hzlog\lfpPCA.mat']);
lf_range = [4,8];
hf_range = [34,72];
twin     = [-5,10]; 
%tPAC_all_session = []; % store PAC mean (PPC1,LPl1,VC1,PPC2,LPl2...)
tPAC_all_session.LPl = [];
tPAC_all_session.PPC = [];
tPAC_all_session.VC = [];
linORlog = 2;
% for irec = 1:numel(recs)
%     display(['rec ' num2str(irec) '\' num2str(numel(recs))])
%     load([recs(irec).folder '\' recs(irec).name]);
%    saveRootPath = recs(irec).folder;
animalCode = '0147';

addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys'));
PreprocessDir = ['E:/FerretData/' animalCode '/Preprocessed/'];
AnalysisDir   = ['E:/FerretData/' animalCode '/Analyzed/'];
BehavDatDir   = ['E:/FerretData/' animalCode '/behav/'];
GroupAnalysisDir = ['E:/FerretData/' animalCode '/GroupAnalysis/PAC/'];
GroupSessionDir = [GroupAnalysisDir 'sessions/'];

fileInfo = dir([PreprocessDir animalCode '_AttentionTask6*']); % detect files to load/convert  '_LateralVideo*'

% loop through each recording
for irec = 1:numel(fileInfo)
    recName = fileInfo(irec).name;
    splitName   = strsplit(recName,'_');
    sessionName = [splitName{2}(14:end) splitName{3}];
    
    %if cluster == 0 && datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180726', 'InputFormat', 'yyyyMMdd'); continue;end

    rootPreprocessDir = [PreprocessDir recName '/'];
    %rootAnalysisDir   = [AnalysisDir recName '/FC' folderSuffix '/']; %eg. FC_median
    % below 3 lines are dealing with multiple record with the same name,
    % use the first record
    BehavNames = dir([BehavDatDir splitName{1} '_Level' splitName{2}(14:end) '_' splitName{3} '_' splitName{4} '*behav.mat']);
    BehavName  = BehavNames(1).name;
    rootBehavDatDir   = [BehavDatDir BehavName];

%     if ~exist(join(rootAnalysisDir),'dir') 
%         mkdir(join(rootAnalysisDir)); fprintf('\nWorking on record %s =============== \n',recName'); end

    % load and process ttl data in ephys
    lfpMat = is_load([rootPreprocessDir 'lfp/lfpMat'],'lfpMat'); % don't feed in denoised data with NaN values
    lfpFs  = is_load([rootPreprocessDir 'lfp/lfpMat'],'lfpFs');
    
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

try
    load([rootPreprocessDir 'eventTimes']);
    evtTime = evtTimes{1};
    lfpMat = is_load([rootPreprocessDir 'lfp\lfpMat.mat'],'lfpMat');
    numRegions = numel(regionNames);
    for iRegion = 1:numRegions
        regionName = regionNames{iRegion};    
        regionlfp = nanmedian(lfpMat(lfp.validChn{iRegion},:),1); % get channel median
        %lfpValid.(regionName)  = lfpMat(lfp.validChn{iRegion},:);
        %% calculate tPAC 
        tempDir = join([GroupSessionDir sessionName '_' regionName '_Init_'],'');
        fprintf(['\nWorking on ' tempDir '\n']); 
        tPAC = PAC_time_resolved(regionlfp,lfpFs,evtTime,twin,tempDir);
        lf_mask = tPAC.lowFreqs>lf_range(1) & tPAC.lowFreqs<lf_range(2);
        hf_mask = tPAC.highFreqs>hf_range(1) & tPAC.highFreqs<hf_range(2);
        tPAC_all_session.(regionName) = [tPAC_all_session.(regionName); reshape(squeeze(nanmean(nanmean(tPAC.plv{1}(hf_mask, lf_mask,:),1),2)),1,[])]; %hf x lf x tvec
        %mean and median looks the same
    end
catch
end
end
save([GroupAnalysisDir 'tPAC_all_session'],'tPAC_all_session');

   

%% plot all sessions

tvec = linspace(twin(1),twin(2),size(tPAC.plv{1},3))' + 1.75;
xLim =[tvec(1),tvec(end)];
regionNames = {'LPl','PPC','VC'};
numRegions  = numel(regionNames);

screensize = get( groot, 'Screensize' );
fig = figure('name','tPAC_all sessions','Position',[10 50 (screensize(3)-100)/2 (screensize(4)-150)/2.5]);
for iRegion = 1:numRegions
    subplot(1,3,iRegion)
    regionName = regionNames{iRegion};
    plot(tvec, tPAC_all_session.(regionName));
    title(regionName); xlabel('Time to init [sec]'); ylabel('theta/gamma PAC');xlim(xLim);
end
savefig(fig, [GroupAnalysisDir 'tPAC_all sessions.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'tPAC_all sessions.png']);


%% plot mean and sem
cm = jet;
ColorSet = [[0,0,1];[1,0,0];[0,0.8,0.2]];

fig = figure('name','tPAC_mean sessions','Position',[10 50 (screensize(3)-100)/2 (screensize(4)-150)/2.5]);
for iRegion = 1:numRegions
    subplot(1,3,iRegion)
    regionName = regionNames{iRegion};
    mean = nanmean(tPAC_all_session.(regionName),1);
    sem  = nanstd(tPAC_all_session.(regionName),[],1)/sqrt(size(tPAC_all_session.(regionName),1));
    shadedErrorBar(tvec, mean, sem, {'color',ColorSet(iRegion,:)}, 0.1);
    title(regionName); xlabel('Time to init [sec]'); ylabel('theta/gamma PAC');
    xlim(xLim);ylim([0.08,0.55]);
end
savefig(fig, [GroupAnalysisDir 'tPAC_mean sessions.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'tPAC_mean sessions.png']);


%% plot median and sem
fig = figure('name','tPAC_mean sessions','Position',[10 50 (screensize(3)-100)/2 (screensize(4)-150)/2.5]);
for iRegion = 1:numRegions
    subplot(1,3,iRegion)
    regionName = regionNames{iRegion};
    median = nanmedian(tPAC_all_session.(regionName),1);
    sem  = nanstd(tPAC_all_session.(regionName),[],1)/sqrt(size(tPAC_all_session.(regionName),1));
    shadedErrorBar(tvec, median, sem,{'color',ColorSet(iRegion,:)}, 0.1);
    title(regionName); xlabel('Time to init [sec]'); ylabel('theta/gamma PAC');
    xlim(xLim);ylim([0.08,0.55]);
end
savefig(fig, [GroupAnalysisDir 'tPAC_mean sessions.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'tPAC_mean sessions.png']);



