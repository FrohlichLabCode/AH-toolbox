clear
clc

animalCode = '0173';

addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys'));
PreprocessDir = ['E:/FerretData/' animalCode '/Preprocessed/'];
AnalysisDir   = ['E:/FerretData/' animalCode '/Analyzed/'];
BehavDatDir   = ['E:/FerretData/' animalCode '/behav/'];
GroupAnalysisDir = ['E:/FerretData/' animalCode '/GroupAnalysis/sessionFC_6b_Delay/'];
GroupSessionDir = [GroupAnalysisDir];
fileInfo = dir([AnalysisDir animalCode '_Level6b*']); % detect files to load/convert  '_LateralVideo*'
alignNames = {'6sDelay'}; %{'Init','Touch'};
alignName = alignNames{1};

condNames   = {'4sDelay', '5sDelay', '6sDelay'};
regionNames = {'FC','LPl','PPC','VC'};
regionPairs = {[1,3],[2,3],[2,4],[3,4]};
numConds    = numel(condNames);
numRegions  = numel(regionNames);
numPairs    = numel(regionPairs);

% generate name pair cells
for iPair = 1:numel(regionPairs)
    regionPairNames{iPair} = [regionNames{regionPairs{iPair}(1)} '-' regionNames{regionPairs{iPair}(2)}];
    regionPair_Names{iPair} = [regionNames{regionPairs{iPair}(1)} '_' regionNames{regionPairs{iPair}(2)}];%{'FC_PPC', 'LPl_PPC', 'LPl_VC', 'PPC_VC'};
    regionPairNamesGC{2*iPair-1} = [regionNames{regionPairs{iPair}(1)} '->' regionNames{regionPairs{iPair}(2)}];%{'LPl->PPC','PPC->LPl';'LPl->VC','VC->LPl'};
    regionPairNamesGC{2*iPair} = [regionNames{regionPairs{iPair}(2)} '->' regionNames{regionPairs{iPair}(1)}];
end

% combine data from all sessions
for irec = 1:numel(fileInfo)
    recName = fileInfo(irec).name;
    splitName   = strsplit(recName,'_');
    for iRegionPair = 1:numel(regionPairNames)
        regionPairName = regionPairNames{iRegionPair};
        regionPair_Name = regionPair_Names{iRegionPair};
        rootAnalysisDir = [AnalysisDir recName '/FC_validChns/' regionPairName '/'];
    for iCond = 1:numConds
        condName = condNames{iCond};
        %load([rootAnalysisDir 'plvAll_' condName '.mat']);
        %try
        if strcmp(regionPairName, 'FC-PPC') 
            Spec.FC(irec,iCond,:,:) = is_load([rootAnalysisDir 'funcCon_median_' alignName '_' condName '.mat'], 'avgXSpec');
            SpecNorm.FC(irec,iCond,:,:) = is_load([rootAnalysisDir 'funcCon_median_' alignName '_' condName '.mat'], 'avgXNormed');
        elseif strcmp(regionPairName, 'LPl-PPC')
            Spec.LPl(irec,iCond,:,:) = is_load([rootAnalysisDir 'funcCon_median_' alignName '_' condName '.mat'], 'avgXSpec');
            SpecNorm.LPl(irec,iCond,:,:) = is_load([rootAnalysisDir 'funcCon_median_' alignName '_' condName '.mat'], 'avgXNormed');
        elseif strcmp(regionPairName, 'PPC-VC')
            Spec.PPC(irec,iCond,:,:) = is_load([rootAnalysisDir 'funcCon_median_' alignName '_' condName '.mat'], 'avgXSpec');
            SpecNorm.PPC(irec,iCond,:,:) = is_load([rootAnalysisDir 'funcCon_median_' alignName '_' condName '.mat'], 'avgXNormed');
            Spec.VC(irec,iCond,:,:) = is_load([rootAnalysisDir 'funcCon_median_' alignName '_' condName '.mat'], 'avgYSpec');
            SpecNorm.VC(irec,iCond,:,:) = is_load([rootAnalysisDir 'funcCon_median_' alignName '_' condName '.mat'], 'avgYNormed');
        end
        
        PLV.(regionPair_Name)(irec,iCond,:,:) = is_load([rootAnalysisDir 'funcCon_median_' alignName '_' condName '.mat'], 'avgPLV');
        Coherence.(regionPair_Name)(irec,iCond,:,:) = abs(is_load([rootAnalysisDir 'funcCon_median_' alignName '_' condName '.mat'], 'avgCoherency'));
        %catch
        %end
        %try
        GC.(regionPair_Name)(irec,iCond,1,:,:) = is_load([rootAnalysisDir 'GC_medain_' alignName '_' condName '.mat'], 'avgGC_XtoY');
        GC.(regionPair_Name)(irec,iCond,2,:,:) = is_load([rootAnalysisDir 'GC_medain_' alignName '_' condName '.mat'], 'avgGC_YtoX');
        if irec==1; GC.tvec = is_load([rootAnalysisDir 'GC_medain_' alignName '_' condName '.mat'],'tvecGC');
            [foi tvec] = is_load([rootAnalysisDir 'specAll_' alignName '_' condName '.mat'],'foi','tvec');end        
        %catch
        %end
    end
    end
end

% shift trial initiation time by 1.75sec
tvec = tvec + 1;
GC.tvec = GC.tvec+1;
[tickLabel,tickLoc] = is_load(['E:\FerretData\0147\GroupAnalysis\Spec_correct\sessions\TrialSpec_6a02_LPl_Init.mat'],'tickLabel','tickLoc');
if ~exist(GroupAnalysisDir,'dir'); mkdir(GroupAnalysisDir);end
save([GroupAnalysisDir 'allSessionFC_byDeley.mat'],'tvec','foi','tickLabel','tickLoc','regionPairNames','regionPair_Names','regionPairNamesGC', 'Spec','SpecNorm','PLV','GC','-v7.3');

%% plot median across sessions
% plot Spec for all regions
fig = figure('name','medianSpec','position', [10 20 320*numConds 270*numRegions]);lw = 2; %x,y,width,height
for iRegion = 1:numRegions %row
    regionName = regionNames{iRegion};
    xLabel = 'Time to init [sec]';
    for iCond = 1:numConds %column
        condName = condNames{iCond};
        subplot(numRegions,numConds,(iRegion-1)*numConds+iCond)
        imagesc(tvec,1:numel(foi),pow2db(squeeze(nanmedian(Spec.(regionName)(:,iCond,:,:),1))));
        xlabel(xLabel); ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([15 45]);
        cl = colorbar('northoutside'); ylabel(cl,[regionName ': ' condName],'FontSize',12)
    end
end
colormap(jet)
savefig(fig, [GroupAnalysisDir 'medianSpec_byDelay.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'medianSpec_byDelay.png']);
% 
fig = figure('name','medianSpecNorm','position', [10 20 320*numConds 270*numRegions]);lw = 2; %x,y,width,height
for iRegion = 1:numRegions %row
    regionName = regionNames{iRegion};
    xLabel = 'Time to init [sec]';
    for iCond = 1:numConds %column
        condName = condNames{iCond};
        subplot(numRegions,numConds,(iRegion-1)*numConds+iCond)
        imagesc(tvec,1:numel(foi),pow2db(squeeze(nanmedian(SpecNorm.(regionName)(:,iCond,:,:),1))));
        xlabel(xLabel); ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-5 5]);
        %if strcmp(regionName, 'FC');caxis([-2 4]);end
        cl = colorbar('northoutside'); ylabel(cl,[regionName ': ' condName],'FontSize',12)
    end
end
colormap(jet)
savefig(fig, [GroupAnalysisDir 'medianSpecNorm_byDelay.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'medianSpecNorm_byDelay.png']);

% plot PLV for 4 pairs
fig = figure('name','medianPLV','position', [10 20 320*numConds 270*numPairs]);lw = 2; %x,y,width,height
for iRegionPair = 1:numPairs
    regionPair_Name = regionPair_Names{iRegionPair};
    xLabel = 'Time to init [sec]';
    % plot PLV
    for iCond = 1:numConds %column
        condName = condNames{iCond};
        subplot(numRegions,numConds,(iRegionPair-1)*numConds+iCond)
        imagesc(tvec,1:numel(foi),squeeze(nanmedian(PLV.(regionPair_Name)(:,iCond,:,:),1)));
        xlabel(xLabel); ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([0.35 0.75]);
        if strcmp(regionPair_Name, 'LPl_PPC');caxis([0.45 0.9]);end
        cl = colorbar('northoutside'); ylabel(cl,[regionPairNames{iRegionPair} ': ' condName],'FontSize',12)
    end
end
colormap(jet)
savefig(fig, [GroupAnalysisDir 'medianPLV_byDelay.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'medianPLV_byDelay.png']);

% plot coherence for 4 pairs
fig = figure('name','medianCoherence','position', [10 20 320*numConds 270*numPairs]);lw = 2; %x,y,width,height
for iRegionPair = 1:numPairs
    regionPair_Name = regionPair_Names{iRegionPair};
    xLabel = 'Time to init [sec]';
    % plot PLV
    for iCond = 1:numConds %column
        condName = condNames{iCond};
        subplot(numPairs,numConds,(iRegionPair-1)*numConds+iCond)
        imagesc(tvec,1:numel(foi),squeeze(nanmedian(Coherence.(regionPair_Name)(:,iCond,:,:),1)));
        xlabel(xLabel); ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([0.05 0.5]);
        %if strcmp(regionPair_Name, 'LPl_PPC');caxis([0.05 0.6]);end
        cl = colorbar('northoutside'); ylabel(cl,[regionPairNames{iRegionPair} ': ' condName],'FontSize',12)
    end
end
colormap(jet)
savefig(fig, [GroupAnalysisDir 'medianCoherence_byDelay.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'medianCoherence_byDelay.png']);


fig = figure('name','medianGC','position', [10 20 320*numConds*2 270*numPairs]);lw = 2; %x,y,width,height
for iRegionPair = 1:numPairs
    regionPair_Name = regionPair_Names{iRegionPair};
    xLabel = 'Time to init [sec]';
    for iCond = 1:numConds %column
        condName = condNames{iCond};
        % X -> Y
        subplot(numPairs,numConds*2,(2*iRegionPair-2)*numConds+iCond)
        imagesc(GC.tvec,1:numel(foi),squeeze(nanmedian(real(GC.(regionPair_Name)(:,iCond,1,:,:)),1)));
        xlabel(xLabel); ylabel('Frequency [Hz]');% title('GC: X to Y')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); ylabel(cl,[regionPairNamesGC{2*iRegionPair-1} ':' condName],'FontSize',12)
        if iRegionPair == 2; caxis([0 0.15]);else caxis([0,0.06]);end 
        ylim([tickLoc(1) tickLoc(end)-10]); % 90+ Hz has saturated values
    
        % Y -> X
        subplot(numPairs,numConds*2,(2*iRegionPair-1)*numConds+iCond)
        imagesc(GC.tvec,1:numel(foi),squeeze(nanmedian(real(GC.(regionPair_Name)(:,iCond,2,:,:)),1)));
        xlabel(xLabel); ylabel('Frequency [Hz]');% title('GC: X to Y')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); ylabel(cl,[regionPairNamesGC{2*iRegionPair} ':' condName],'FontSize',12)
        if iRegionPair == 2; caxis([0 0.15]);else caxis([0,0.06]);end
        ylim([tickLoc(1) tickLoc(end)-10]); % 90+ Hz has saturated values   
    end
end
colormap(jet)
savefig(fig, [GroupAnalysisDir 'medianGC_byDelay.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'medianGC_byDelay.png']);

close all
        