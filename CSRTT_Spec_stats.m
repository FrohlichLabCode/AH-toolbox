% Spec mean across trials and sessions

clear
tic

cluster = 0;
skipRec = 1;
linORlog = 2; %freqs of interest: 1=linear 2=log
MedianorPCA = 3; 
animalCode = '0147';
regionNames = {'LPl','PPC','VC'};
numRegions  = numel(regionNames);
condNames   = {'Init'};%,'Touch'};
numConds    = numel(condNames);
twins = {[-5 10];[-4 5]};
foiNames = {'Theta','Alpha','Gamma'};
foiWins = {[4,8],[10,14],[32,70]}; % from theta-gamma coupling plot
numFreqs  = numel(foiNames);

%cm = jet;
%ColorSet = cm([60,50,30],:);


ColorSet = [[0,0,205];[138,43,226];[238,130,238]]/256; %medianblue, blueviolet, violet, from https://www.rapidtables.com/web/color/purple-color.html
%load('E:\FerretData\0147\Preprocessed\0147_AttentionTask6a_04_20170927\eventTimes'); 
[baseTwins,foi,fois,tickLabel,tickLoc,t] = is_load(['E:\FerretData\0147\GroupAnalysis\Spec_correct\sessions\TrialSpec_6a02_LPl_Init.mat'],'baseTwins','foi','fois','tickLabel','tickLoc','tvec');
tvec{1} = t + 1.75; % offset of beam pass
%tvec{2} = is_load(['E:\FerretData\0147\GroupAnalysis\Spec_missed\sessions\TrialSpec_6a02_LPl_Touch.mat'],'tvec');
for iFreq = 1:numFreqs
    foiMask(iFreq,:) = foi>= foiWins{iFreq}(1) & foi<= foiWins{iFreq}(2);
end


if cluster == 0 %linux use '/', windows matlab can use both '/' and '\'
    PreprocessDir = ['E:/FerretData/' animalCode '/Preprocessed/'];
    AnalysisDir   = ['E:/FerretData/' animalCode '/Analyzed/'];
    BehavDatDir   = ['E:/FerretData/' animalCode '/behav/'];
    GroupAnalysisDir = ['E:/FerretData/' animalCode '/GroupAnalysis/Spec_premature/'];
elseif cluster == 1
    addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
    PreprocessDir = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Preprocessed/']; % CHANGE FOR KILLDEVIL VS LONGLEAF
    AnalysisDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Analyzed/'];
    BehavDatDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/behav/'];
    GroupAnalysisDir = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/GroupAnalysis/Spec/'];

    %code for initialising parallel computing
%     numCore = 24; % USR DEFINE
%     myPool = parpool('local',numCore,'SpmdEnabled',false);  
end

trialStartInd = ones(numConds,numRegions);
trialEndInd = zeros(numConds,numRegions);


% Pool all trials of one region together
for iCond = 1:numConds
    condName = condNames{iCond};
    for iRegion = 1:numRegions
        regionName = regionNames{iRegion};
        %trialPowerAll.(regionName) = []; %
        fileInfo = dir([GroupAnalysisDir 'sessions/TrialSpec_*' regionName '_' condName '.mat']); % detect files to load/convert  '_LateralVideo*'
        for iSession = 1:numel(fileInfo)
            fileName = fileInfo(iSession).name;
            fprintf(['Processing session ' fileName '\n']);
            spec = is_load([GroupAnalysisDir 'sessions/' fileName],'TrialSpec');
            if length(size(spec.(regionName)))== 3; numTrials = size(spec.(regionName),1);
                elseif length(size(spec.(regionName))) ==2; numTrials = 1;end
            if iCond == 1                 
                trialEndInd(iCond,iRegion) = trialStartInd(iCond,iRegion)+numTrials-1;
                trialPowerAll_Init.(regionName)(trialStartInd(iCond,iRegion):trialEndInd(iCond,iRegion),:,:) = spec.(regionName);
                trialPowerAllnormed_Init.(regionName)(trialStartInd(iCond,iRegion):trialEndInd(iCond,iRegion),:,:) = spec.([regionName '_normed']);
                trialStartInd(iCond,iRegion) = trialEndInd(iCond,iRegion)+1;
            elseif iCond == 2
                trialEndInd(iCond,iRegion) = trialStartInd(iCond,iRegion)+numTrials-1;
                trialPowerAll_Touch.(regionName)(trialStartInd(iCond,iRegion):trialEndInd(iCond,iRegion),:,:) = spec.(regionName);
                trialPowerAllnormed_Touch.(regionName)(trialStartInd(iCond,iRegion):trialEndInd(iCond,iRegion),:,:) = spec.([regionName '_normed']);
                trialStartInd(iCond,iRegion) = trialEndInd(iCond,iRegion)+1;
            end            
        end
    end
end

%save([GroupAnalysisDir 'trialPool_PowerAll.mat'],'trialPowerAll_Init','trialPowerAllnormed_Init','trialPowerAll_Touch','trialPowerAllnormed_Touch', '-v7.3');
save([GroupAnalysisDir 'trialPool_PowerAll.mat'],'trialPowerAll_Init','trialPowerAllnormed_Init', '-v7.3');

%% first visualize all trial average spectrogram for init and touch
% median
for iCond = 1:numConds % stimulation frequency
    condName = condNames{iCond};
    fig = figure('name',[condName ': trialPool_MedianSpec'],'position', [10   20   320*numRegions   270*2]);lw = 2; %x,y,width,height
    xLabel=(['Time to ' condName ' [sec]']);
    for iRegion = 1:numRegions
        regionName = regionNames{iRegion};
        subplot(2,numRegions,iRegion)
        if iCond == 1
        mat2plot = pow2db(squeeze(nanmedian(trialPowerAll_Init.(regionName),1)));
        elseif iCond == 2
        mat2plot = pow2db(squeeze(nanmedian(trialPowerAll_Touch.(regionName),1)));
        end
        imagesc(tvec{iCond},1:numel(foi),mat2plot);
        xlabel(xLabel); ylabel('Frequency [Hz]'); % title('Signal X power')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        %caxis([-2 2]);
        cl = colorbar('northoutside'); 
        ylabel(cl,['Power [dB]: ' regionName],'FontSize',12)
    
        subplot(2,numRegions,iRegion+numRegions)
        if iCond == 1
        mat2plot = pow2db(squeeze(nanmedian(trialPowerAllnormed_Init.(regionName),1)));
        elseif iCond == 2
        mat2plot = pow2db(squeeze(nanmedian(trialPowerAllnormed_Touch.(regionName),1)));
        end
        imagesc(tvec{iCond},1:numel(foi),mat2plot);
        xlabel(xLabel); ylabel('Frequency [Hz]'); % title('Signal X power')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        %caxis([-2 2]); 
        cl = colorbar('northoutside'); ylabel(cl,['Normed power [dB]: ' regionName],'FontSize',12)
    end
    colormap(jet)
    savefig(fig, [GroupAnalysisDir 'trialPool_MedianPower_' condName '.fig'],'compact');
    saveas(fig, [GroupAnalysisDir 'trialPool_MedianPower_' condName '.png']);
end

%% same for mean -- looks bad, contaminated by artifacts, so just use median
% for iCond = 1:numConds % stimulation frequency
%     condName = condNames{iCond};
%     fig = figure('name',[condName ': trialPool_MeanSpec'],'position', [10   20   320*numRegions   270*2]);lw = 2; %x,y,width,height
%     for iRegion = 1:numRegions
%         regionName = regionNames{iRegion};
%         subplot(2,numRegions,iRegion)
%         if iCond == 1
%         mat2plot = pow2db(squeeze(nanmean(trialPowerAll_Init.(regionName),1)));
%         elseif iCond == 2
%         mat2plot = pow2db(squeeze(nanmean(trialPowerAll_Touch.(regionName),1)));
%         end
%         imagesc(tvec{iCond},1:numel(foi),mat2plot);
%         xlabel(xLabel); ylabel('Frequency [Hz]'); % title('Signal X power')
%         ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
%         %caxis([-2 2]); 
%         cl = colorbar('northoutside'); ylabel(cl,['Power [dB]: ' regionName],'FontSize',12)
%         
%     
%         subplot(2,numRegions,iRegion+numRegions)
%         if iCond == 1
%         mat2plot = pow2db(squeeze(nanmean(trialPowerAllnormed_Init.(regionName),1)));
%         elseif iCond == 2
%         mat2plot = pow2db(squeeze(nanmean(trialPowerAllnormed_Touch.(regionName),1)));
%         end
%         imagesc(tvec{iCond},1:numel(foi),mat2plot);
%         xlabel(xLabel); ylabel('Frequency [Hz]'); % title('Signal X power')
%         ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
%         %caxis([-2 2]); 
%         cl = colorbar('northoutside'); ylabel(cl,['Normed power [dB]: ' regionName],'FontSize',12)
%     end
%     colormap(jet)
%     savefig(fig, [GroupAnalysisDir 'trialPool_MeanSpec_' condName '.fig'],'compact');
%     saveas(fig, [GroupAnalysisDir 'trialPool_MeanSpec_' condName '.png']);
% end



%% Calculate foi average
for iRegion = 1:numRegions
    regionName = regionNames{iRegion};
    for iFreq = 1:numFreqs
        foiName = foiNames{iFreq};
        for iCond = 1:numConds
            condName = condNames{iCond};
            if strcmp(condName,'Init')
            foiPowerAll_Init.(regionName)(iFreq,:,:) = squeeze(nanmean(trialPowerAll_Init.(regionName)(2:end,foiMask(iFreq,:),:),2));
            foiPowerAllnormed_Init.(regionName)(iFreq,:,:) = squeeze(nanmean(trialPowerAllnormed_Init.(regionName)(2:end,foiMask(iFreq,:),:),2));
            elseif strcmp(condName,'Touch')
            foiPowerAll_Touch.(regionName)(iFreq,:,:) = squeeze(nanmean(trialPowerAll_Touch.(regionName)(2:end,foiMask(iFreq,:),:),2));
            foiPowerAllnormed_Touch.(regionName)(iFreq,:,:) = squeeze(nanmean(trialPowerAllnormed_Touch.(regionName)(2:end,foiMask(iFreq,:),:),2));
    
            end
        end
    end
end
%save([GroupAnalysisDir 'trialPool_foiPowerAll.mat'],'foiPowerAll_Init','foiPowerAllnormed_Init','foiPowerAll_Touch','foiPowerAllnormed_Touch','-v7.3');
save([GroupAnalysisDir 'trialPool_foiPowerAll.mat'],'foiPowerAll_Init','foiPowerAllnormed_Init','-v7.3');

% iFreq x iTrial x tvec

%% Plot time-resolved foi average Spec
for iCond = 1%:numConds
    condName = condNames{iCond};
    fig = figure('name',[condName ': trialPool_foiSpec'],'position', [10   20   320*numRegions   270*2]);lw = 2; %x,y,width,height
    xLabel=(['Time to ' condName ' [sec]']); xLim = {[-2,10];[-4,5]};
    for iRegion = 1:numRegions
        regionName = regionNames{iRegion};
        subplot(2,numRegions,iRegion)
        for iFreq = 1:numFreqs
            foiName = foiNames{iFreq};
            if iCond == 1
            median = squeeze(pow2db(nanmedian(foiPowerAll_Init.(regionName)(iFreq,:,:),2)))';
            sem    = squeeze(pow2db(nanstd(foiPowerAll_Init.(regionName)(iFreq,:,:),[],2)))'/sqrt(size(foiPowerAll_Init.(regionName),2));            
            elseif iCond == 2
            median = squeeze(pow2db(nanmedian(foiPowerAll_Touch.(regionName)(iFreq,:,:),2)))';
            sem    = squeeze(pow2db(nanstd(foiPowerAll_Touch.(regionName)(iFreq,:,:),[],2)))'/sqrt(size(foiPowerAll_Touch.(regionName),2));
            end
            H(iFreq) = shadedErrorBar(tvec{iCond}, median, sem,{'color',ColorSet(iFreq,:),'LineWidth',4}, 0.1); hold on; %last is transparency level            
%            H(iFreq) = shadedErrorBar(tvec, median, sem,{'color',ColorSet(iFreq,:),'LineWidth',4}, 0.1); hold on; %last is transparency level            
            
        end
        title(regionName);xlim(xLim{iCond});
        if iRegion == 1; legend([H(1).mainLine H(2).mainLine H(3).mainLine], 'theta','alpha','gamma'); clear H; end
        xlabel(xLabel); ylabel('Power [dB]'); % title('Signal X power')
        xticks([-2,0,5,10]);
        
        % plot normalized version
        subplot(2,numRegions,iRegion+numRegions)
        for iFreq = 1:numFreqs
            foiName = foiNames{iFreq};
            if iCond == 1
            median = squeeze(pow2db(nanmedian(foiPowerAllnormed_Init.(regionName)(iFreq,:,:),2)));
            sem    = squeeze(pow2db(nanstd(foiPowerAllnormed_Init.(regionName)(iFreq,:,:),[],2)))/sqrt(size(foiPowerAllnormed_Init.(regionName),2));
            elseif iCond == 2
            median = squeeze(pow2db(nanmedian(foiPowerAllnormed_Touch.(regionName)(iFreq,:,:),2)));
            sem    = squeeze(pow2db(nanstd(foiPowerAllnormed_Touch.(regionName)(iFreq,:,:),[],2)))/sqrt(size(foiPowerAllnormed_Touch.(regionName),2));                        
            end
            shadedErrorBar(tvec{iCond}, median, sem,{'color',ColorSet(iFreq,:)}, 0.1); hold on; %last is transparency level            
        end
        xlim(xLim{iCond});xlabel(xLabel); ylabel('BLNormed power [dB]'); % title('Signal X power')
        xticks([-2,0,5,10]);
    end
    fig.Color = 'white';
    savefig(fig, [GroupAnalysisDir 'trialPool_foiSpec_' condName '.fig'],'compact');
    saveas(fig, [GroupAnalysisDir 'trialPool_foiSpec_' condName '.png']);
end

