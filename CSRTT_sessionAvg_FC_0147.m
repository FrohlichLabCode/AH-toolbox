animalCode = '0147';

addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys'));
PreprocessDir = ['E:/FerretData/' animalCode '/Preprocessed/'];
AnalysisDir   = ['E:/FerretData/' animalCode '/Analyzed/'];
BehavDatDir   = ['E:/FerretData/' animalCode '/behav/'];
GroupAnalysisDir = ['E:/FerretData/' animalCode '/GroupAnalysis/FC/'];
GroupSessionDir = [GroupAnalysisDir 'sessions/'];
fileInfo = dir([AnalysisDir animalCode '_AttentionTask6b*']); % detect files to load/convert  '_LateralVideo*'

condNames = {'Init','Touch'};

for irec = 1:numel(fileInfo)
    recName = fileInfo(irec).name;
    splitName   = strsplit(recName,'_');
    rootAnalysisDir = [AnalysisDir recName '/FC_validChns_150f/LPl-PPC/'];
    for iCond = 1:numConds
        condName = condNames{iCond};
        %load([rootAnalysisDir 'plvAll_' condName '.mat']);
        try
        PLV.LPl_PPC(irec,:,:) = is_load([rootAnalysisDir 'funcCon_avg_' condName '.mat'], 'avgPLV');
        catch
        end
        try
        GC.LPl_PPC(irec,1,:,:) = is_load([rootAnalysisDir 'GCAll_' condName '.mat'], 'avgGC_XtoY');
        GC.LPl_PPC(irec,2,:,:) = is_load([rootAnalysisDir 'GCAll_' condName '.mat'], 'avgGC_YtoX');
        if irec==1; GC.tvec(iCond,:) = is_load([rootAnalysisDir 'GCAll_' condName '.mat'],'tvecGC');
            [PLV.foi PLV.tvec(iCond,:)] = is_load([rootAnalysisDir 'specAll_' condName '.mat'],'foi','tvec');end        
        catch
        end
   end
end
% shift trial initiation time by 1.75sec
PLV.tvec(1,:) = PLV.tvec(1,:)+1.75;
GC.tvec(1,:)  = GC.tvec(1,:)+1.75;

fileInfo = dir([AnalysisDir animalCode '_AttentionTask6a*']); % detect files to load/convert  '_LateralVideo*'
for irec = 1:numel(fileInfo)
    recName = fileInfo(irec).name;
    splitName   = strsplit(recName,'_');
    rootAnalysisDir = [AnalysisDir recName '/FC_validChns_150f/LPl-VC/'];
    for iCond = 1:numConds
        condName = condNames{iCond};
        %load([rootAnalysisDir 'plvAll_' condName '.mat']);
        try
        PLV.LPl_VC(irec,:,:) = is_load([rootAnalysisDir 'funcCon_avg_' condName '.mat'], 'avgPLV');
        catch
        end
        try
        GC.LPl_VC(irec,1,:,:) = is_load([rootAnalysisDir 'GCAll_' condName '.mat'], 'avgGC_XtoY');
        GC.LPl_VC(irec,2,:,:) = is_load([rootAnalysisDir 'GCAll_' condName '.mat'], 'avgGC_YtoX');
        catch
        end
   end
end

save([GroupAnalysisDir 'allSessionFC.mat'],'PLV','GC', '-v7.3');

regionPairs = {'LPl-PPC','LPl-VC'};
region_Pairs = {'LPl_PPC','LPl_VC'};
region_PairsGC = {'LPl->PPC','PPC->LPl';'LPl->VC','VC->LPl'};
[tickLabel,tickLoc] = is_load(['E:\FerretData\0147\GroupAnalysis\Spec\sessions\TrialSpec_6a02_LPl_Init.mat'],'tickLabel','tickLoc');

% plot FC for 2 pairs
PLV.LPl_PPC = PLV.LPl_PPC(1:4,:,:);
GC.LPl_PPC  = GC.LPl_PPC(1:4,:,:,:);
fig = figure('name','medianFC_PLV_GC','position', [10 20 320*3 270*2]);lw = 2; %x,y,width,height
for iRegionPair = 1:numel(regionPairs)
    region_Pair = region_Pairs{iRegionPair};
    xLabel = 'Time to init [sec]';
    % plot PLV
    subplot(2,3,iRegionPair*3-2)
    imagesc(PLV.tvec,1:numel(foi),squeeze(nanmedian(PLV.(region_Pair),1)));    % only first 4 sessions look correct   
    xlabel(xLabel); ylabel('Frequency [Hz]');% title('PLV')
    ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    %caxis([0.1 0.7]);
    cl = colorbar('northoutside'); ylabel(cl,['PLV: ' regionPairs{iRegionPair}],'FontSize',12)
    
    % plot GC X to Y
    subplot(2,3,iRegionPair*3-1)
    imagesc(GC.tvec,1:numel(foi),squeeze(nanmedian(real(GC.(region_Pair)(:,1,:,:)),1)));
    xlabel(xLabel); ylabel('Frequency [Hz]');% title('GC: X to Y')
    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    %cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: X to Y','FontSize',15)
    cl = colorbar('northoutside'); ylabel(cl,['GC: ' region_PairsGC{iRegionPair,1}],'FontSize',12)
    if iRegionPair == 1; caxis([0 1]);else caxis([0,0.3]);end 
    ylim([tickLoc(1) tickLoc(end)-10]); % 90+ Hz has saturated values
    
    % plot GC Y to X
    subplot(2,3,iRegionPair*3)
    imagesc(GC.tvec,1:numel(foi),squeeze(nanmedian(real(GC.(region_Pair)(:,2,:,:)),1)));
    xlabel(xLabel); ylabel('Frequency [Hz]');% title('GC: X to Y')
    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    %cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: X to Y','FontSize',15)
    cl = colorbar('northoutside'); ylabel(cl,['GC: ' region_PairsGC{iRegionPair,2}],'FontSize',12)
    if iRegionPair == 1; caxis([0 1]);else caxis([0,0.3]);end 
    ylim([tickLoc(1) tickLoc(end)-10]); % 90+ Hz has saturated values
    
end
colormap(jet)

savefig(fig, [GroupAnalysisDir 'medianFuncCon_PLV_GC.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'medianFuncCon_PLV_GC.png']);


%% plot mean instead, looks smoother
fig = figure('name','meanFC_PLV_GC','position', [10 20 320*3 270*2]);lw = 2; %x,y,width,height
for iRegionPair = 1:numel(regionPairs)
    region_Pair = region_Pairs{iRegionPair};
    xLabel = 'Time to init [sec]';
    % plot PLV
    subplot(2,3,iRegionPair*3-2)
    imagesc(PLV.tvec,1:numel(foi),squeeze(nanmean(PLV.(region_Pair),1)));    % only first 4 sessions look correct   
    xlabel(xLabel); ylabel('Frequency [Hz]');% title('PLV')
    ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    %if iRegionPair == 1; caxis([0.1 0.9]);else caxis('auto');end
    cl = colorbar('northoutside'); ylabel(cl,['PLV: ' regionPairs{iRegionPair}],'FontSize',12)
    set(gca,'XTick',[0,5,10]);
    
    % plot GC X to Y
    subplot(2,3,iRegionPair*3-1)
    imagesc(GC.tvec,1:numel(foi),squeeze(nanmean(real(GC.(region_Pair)(:,1,:,:)),1)));
    xlabel(xLabel); ylabel('Frequency [Hz]');% title('GC: X to Y')
    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    %cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: X to Y','FontSize',15)
    cl = colorbar('northoutside'); ylabel(cl,['GC: ' region_PairsGC{iRegionPair,1}],'FontSize',12)
    if iRegionPair == 1; caxis([0 1]);else caxis([0,0.3]);end 
    ylim([tickLoc(1) tickLoc(end)-10]); % 90+ Hz has saturated values
    set(gca,'XTick',[0,5,10]);
    
    % plot GC Y to X
    subplot(2,3,iRegionPair*3)
    imagesc(GC.tvec,1:numel(foi),squeeze(nanmean(real(GC.(region_Pair)(:,2,:,:)),1)));
    xlabel(xLabel); ylabel('Frequency [Hz]');% title('GC: X to Y')
    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    %cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: X to Y','FontSize',15)
    cl = colorbar('northoutside'); ylabel(cl,['GC: ' region_PairsGC{iRegionPair,2}],'FontSize',12)
    if iRegionPair == 1; caxis([0 1]);else caxis([0,0.3]);end 
    ylim([tickLoc(1) tickLoc(end)-10]); % 90+ Hz has saturated values
    set(gca,'XTick',[0,5,10]);
end
colormap(jet)

savefig(fig, [GroupAnalysisDir 'meanFuncCon_PLV_GC.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'meanFuncCon_PLV_GC.png']);



%% plot PLV spectra before, during and after touch
timeWin = [[-2,0]; [0,4]; [9,10]];
timeWinNames = {'Before initiation','During delay','After touch'};
for itimeWin = 1:size(timeWin,1)
    tvecMask{itimeWin} = PLV.tvec>= timeWin(itimeWin,1) & PLV.tvec <= timeWin(itimeWin,2);
end

fig = figure('name','meanPLV_3timeWin','position', [10 20 320 270*2]);lw = 4; %x,y,width,height
[tickLoc, tickLabel] = getTickLabel(lowFreq, highFreq, numFreqs, 2);

for iRegionPair = 1:numel(regionPairs)
    region_Pair = region_Pairs{iRegionPair};
    regionPair  = regionPairs{iRegionPair};
    for itimeWin = 1:size(timeWin,1)
        timeWinName = timeWinNames{itimeWin};

        subplot(2,1,iRegionPair)
        spectra = squeeze(nanmean(PLV.(region_Pair)(:,:,tvecMask{itimeWin}),3));
        mean = squeeze(nanmean(spectra,1));
        sem  = squeeze(nanstd(spectra,[],1));
        semilogx(foi, mean,'linewidth',lw);hold on; 
        %try to change tick label to 2,4,8... but unsuccessful
%         set(gca, 'XTickLabel',[]);                      %# suppress current x-labels
%         xt = get(gca, 'XTick');yl = get(gca, 'YLim');
%         hTxt = text(xt, yl(ones(size(xt))), tickLabel, ...   %# create text at same locations
%         'Interpreter','tex', ...                   %# specify tex interpreter
%         'VerticalAlignment','top', ...             %# v-align to be underneath
%         'HorizontalAlignment','center');
        title(regionPair);xlim([2,96]);
        %xlim([tickLoc(1),tickLoc(end)-10]);
        if itimeWin == 3; set(gca,'XDir','reverse');camroll(-90);end
    end
    legend(timeWinNames{:});
end

savefig(fig, [GroupAnalysisDir 'meanPLV_3timeWin.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'meanPLV_3timeWin.png']);



close all
        