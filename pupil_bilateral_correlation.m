% an example session
rootPreprocessDir = 'D:\FerretData\0169\Preprocessed\0169_LateralVideo_033_20180808\';
rootAnalysisDir   = 'D:\FerretData\0169\Analyzed\0169_LateralVideo_032_20180808\Pupil\';
excludeWin = [-0.1 0.6]; % in seconds

[evtTimes, condNames, condID] = is_load([rootPreprocessDir 'eventTimes'], 'evtTimes', 'condNames', 'condID');
lfpFs = is_load([rootPreprocessDir 'lfp/lfpMat'],'lfpFs');
load([rootPreprocessDir ,'pupil']); % directly load saved pupil data
eyeData = [cazi;cele;cpup]; % flip so that left eye is first, right eye is second
% adc channel names
numPul = size(eyeData,1);
eyeChn = {'leftX','rightX','leftY','rightY','leftD','rightD'};

% get rid of stark outliers
elimThresh = [nanmedian(eyeData,2)-2.5*nanstd(eyeData,0,2) nanmedian(eyeData,2)+2.5*nanstd(eyeData,0,2)]; % eliminate outliers for this analysis
eyeData(eyeData<elimThresh(:,1)) = nan; % get rid out stark outliers
eyeData(eyeData>elimThresh(:,2)) = nan; % get rid out stark outliers

% collapse data for diff condition

for iCond = 1:numel(condID)
    evtMask{iCond} = logical(zeros(numel(condID),size(eyeData,2)));
    for ievt = 1:numel(evtTimes{iCond})
        evtStart = round((evtTimes{iCond}(ievt) + excludeWin(2))*lfpFs);
        evtEnd   = round((evtTimes{iCond}(ievt) + excludeWin(2) + 5)*lfpFs);
        evtMask{iCond}(:,evtStart:evtEnd) = true;
    end
end

screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-150)/2 (screensize(4)-150)]);

for i = 1:3
    subplot(3,3,i)
    vecA = eyeData(2*i-1,:); %60*rawFs:2*60*rawFs);
    vecB = eyeData(2*i,:); %60*rawFs:2*60*rawFs);
    keepMask = ~isnan(vecA) & ~isnan(vecB);
    scatter(vecA(keepMask),vecB(keepMask)); %
    R = corrcoef(vecA(keepMask),vecB(keepMask));
    lsline;
    title(['All trials: Corr coef = ' num2str(R(1,2))]);
    xlabel(eyeChn{2*i-1});ylabel(eyeChn{2*i});
    
    subplot(3,3,i+3)
    vecA = eyeData(2*i-1,evtMask{1}(1,:)); % left video trials
    vecB = eyeData(2*i,evtMask{1}(2,:)); 
    keepMask = ~isnan(vecA) & ~isnan(vecB);
    scatter(vecA(keepMask),vecB(keepMask)); %
    R = corrcoef(vecA(keepMask),vecB(keepMask));
    lsline;
    title(['Left video trials: Corr coef = '  num2str(R(1,2))]);
    xlabel(eyeChn{2*i-1});ylabel(eyeChn{2*i});
    
    subplot(3,3,i+6)
    vecA = eyeData(2*i-1,evtMask{2}(1,:)); %60*rawFs:2*60*rawFs);
    vecB = eyeData(2*i,evtMask{2}(2,:)); %60*rawFs:2*60*rawFs);
    keepMask = ~isnan(vecA) & ~isnan(vecB);
    scatter(vecA(keepMask),vecB(keepMask)); %
    R = corrcoef(vecA(keepMask),vecB(keepMask));
    lsline;
    title(['Right video trials: Corr coef = '  num2str(R(1,2))]);
    xlabel(eyeChn{2*i-1});ylabel(eyeChn{2*i});    
    
end

savefig(fig, [rootAnalysisDir 'bilateral eye correlation.fig'],'compact');
saveas(fig, [rootAnalysisDir 'bilateral eye correlation.png']);
