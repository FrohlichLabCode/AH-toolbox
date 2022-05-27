clear
% This script will analyze an entire session (resting state) and plot
% power spectra

addpath(genpath('E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
skipRec = 1;
% set linear or log plotting scale
numFreqs = 150;
[foi, tickLoc, tickLabel] = getFoiLabel(2, 128, numFreqs, 2); % (lowFreq, highFreq, numFreqs, linORlog)

animalCodes = {'0182','0185','0183'};

for iAnimal = 3%1:numel(animalCodes)
    animalCode = animalCodes{iAnimal};
    
%     PreprocessDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/Preprocessed/'];
%     AnalysisDir   = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/Analyzed/'];
    PreprocessDir = ['Z:/Individual/Angel/FerretData/' animalCode '/Preprocessed/'];
    GroupAnalysisDir   = ['Z:/Individual/Angel/FerretData/' animalCode '/GroupAnalysis_Resting/'];
    %BehavDatDir   = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/behav/'];
    fileInfo      = dir([PreprocessDir animalCode '_Resting*']);  %_Level6* detect files to load/convert  '_LateralVideo*'
    
    if exist([GroupAnalysisDir 'spectra*']) % skip already analyzed records
        fprintf('Animal %s already analyzed \n', animalCode); 
        if skipRec == 1; continue; end
    else
        fprintf('Analyzing animal %s \n',animalCode); 
    end
    
    % region info
    if ismember(animalCode, {'0182','0183','0185'})
        regionNames = {'PPC'};
        numRegion   = numel(regionNames);
        allChn      = {[1:16]};
    end
    numSess = numel(fileInfo);
    if exist([GroupAnalysisDir 'GroupSpectra_' num2str(numSess) 'Sess.mat'])
        load([GroupAnalysisDir 'GroupSpectra_' num2str(numSess) 'Sess.mat']);
    else
    % Initialize empty matrix
    for iRegion = 1:numRegion
        numChn = numel(allChn{iRegion});
        regionName = regionNames{iRegion};
        spectra.(regionName) = NaN(numel(fileInfo),numChn,numFreqs);
    end
    % loop through each recording
    lastsize = 0;
            
    for irec = 1:numel(fileInfo)
        recName = fileInfo(irec).name;   %recName = '0168_Opto_010_20180713';
        fprintf(repmat('\b', 1, lastsize));
        lastsize = fprintf(['Processing ' recName]);

        splitName = strsplit(recName,'_');
        %if datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180823', 'InputFormat', 'yyyyMMdd'); continue;end

        rootPreprocessDir = [PreprocessDir recName '/'];
        %rootAnalysisDir   = [AnalysisDir recName '/LFP/'];
        %rootBehavDatDir   = [BehavDatDir recName '/'];  
        
        [lfpMat,lfpFs] = is_load([rootPreprocessDir 'lfp/lfpMat'],'lfpMat','lfpFs'); % including lfpFs=1000    
        if irec == 1 % only run once
            morWav = sub_makeWavelet(foi,lfpFs);
        end
        %% extract snippits of lfp
        xRange = [10*lfpFs, size(lfpMat,2)-10*lfpFs]; % exclude beginning and end 10sec. in sec, window around 2nd trial Init twins [-5,11]
        xMask = false(1,size(lfpMat,2));
        xMask(xRange(1):xRange(2)) = true;
        tRange = round(xRange(1)*lfpFs):round(xRange(end)*lfpFs); %column index
        for iRegion = 1:numRegion
            numChn = numel(allChn{iRegion});
            regionName = regionNames{iRegion};
            
            for iChn = 1:numChn
                xser = lfpMat(allChn{iRegion}(iChn),xMask);
                xsig = sub_rejectNoise(xser,lfpFs,2,1);
                xspec = nan(numFreqs,numel(xsig));
                dispstat('','init'); % One time only initialization
                for f = 1:numFreqs
                    xspec(f,:) = conv(xsig,morWav{f},'same');
                end
                xmat = nanmean(xspec, 2); % get spectra    
                xpow = abs(xmat).^2;
                spectra.(regionName)(irec,iChn,:) = xpow; 
                spectraNorm.(regionName)(irec,iChn,:) = xpow./sum(xpow);
                spectraNormdb.(regionName)(irec,iChn,:) = pow2db(xpow./sum(xpow));
            end
        end
    end
    AH_mkdir(GroupAnalysisDir);
    numSess = size(spectra.(regionName),1);
    save([GroupAnalysisDir 'GroupSpectra_' num2str(numSess) 'Sess.mat'],'spectra','spectraNorm','spectraNormdb','foi','xRange','tRange','fileInfo', '-v7.3');
    end
    
    %% Plot entire recording spectra, 1 figure per region
    keepSess = 1:numSess;
    if strcmp(animalCode, '0182') % visually exclude sessions
        keepSess = [3:4,10:25]; 
    elseif strcmp(animalCode, '0185')
        keepSess = [2:9,11:17,19:27];
    elseif strcmp(animalCode, '0183')
        keepSess = [1:numSess];
    end
    
    for iRegion = 1:numRegion
        numChn = numel(allChn{iRegion});
        regionName = regionNames{iRegion};
        % select sessions
        spectraNormdb.(regionName) = spectraNormdb.(regionName)(keepSess,:,:);
        
        numSess = size(spectraNormdb.(regionName),1);
        mnChn = squeeze(nanmean(spectraNormdb.(regionName),2));
        mdChn = squeeze(nanmedian(spectraNormdb.(regionName),2));
        
        saveName = ['GroupSpectra_' regionName '_' num2str(numSess) 'Sess'];
        xLim = [50,76]; % 26=4Hz; 50=8Hz; 76=16Hz; 100=32Hz
        yLabel = 'PowerNorm [dB]';
        fig = AH_figure(3,2,saveName);
        subplot(321)
        plot(1:numel(foi), mnChn)
        set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
        %set(gca,'YScale','log')
        title([animalCode ' n=' num2str(numSess) ' allSess mnChn']);
        xlabel('Frequency [Hz]'); xlim(xLim);
        ylabel(yLabel);               
        
        subplot(322)
        plot(1:numel(foi), mdChn); 
        set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
        %set(gca,'YScale','log')
        title(['n=' num2str(numSess) ' allSess mdChn']);
        xlabel('Frequency [Hz]'); xlim(xLim);
        ylabel(yLabel);
        
        subplot(323)
        plot(1:numel(foi), mnChn+pow2db(foi))
        set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
        %set(gca,'YScale','log')
        title(['n=' num2str(numSess) ' allSess mnChn']);
        xlabel('Frequency [Hz]'); xlim(xLim);
        ylabel('PowerNorm*foi [dB]');        
        
        subplot(324)
        plot(1:numel(foi), mdChn+pow2db(foi)); 
        set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
        %set(gca,'YScale','log')
        title(['n=' num2str(numSess) ' allSess mdChn']);
        xlabel('Frequency [Hz]'); xlim(xLim);
        ylabel('PowerNorm*foi [dB]');
        
        clear alphaIDText
        if strcmp(animalCode,'0182')
            alphaID = [62, 66, 69];
        elseif strcmp(animalCode,'0185')
            alphaID = [56, 64, 67, 73];
        elseif strcmp(animalCode,'0183')
            alphaID = [60, 66];    
        end
        for i = 1:numel(alphaID)
            alphaIDText{i} = [num2str(round(foi(alphaID(i)),2)) 'Hz'];
        end
        subplot(325)
        shadedErrorBar(1:numel(foi),mnChn,{@mean,@std},{'k-'}); 
        set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
        %set(gca,'YScale','log')
        title(['n=' num2str(numSess) ' mnSess mnChn + std']);
        xlabel('Frequency [Hz]'); xlim(xLim);
        ylabel(yLabel);
        vline(alphaID);
        text(alphaID, mean(mnChn(:,alphaID),1)+1, alphaIDText)
        
        subplot(326)
        shadedErrorBar(1:numel(foi),mdChn,{@mean,@std},{'k-'}); 
        set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
        %set(gca,'YScale','log')
        title(['n=' num2str(numSess) ' mnSess mdChn + std']);
        xlabel('Frequency [Hz]'); xlim(xLim);
        ylabel(yLabel);
        vline(alphaID);
        text(alphaID, mean(mdChn(:,alphaID),1)+1, alphaIDText)
        
        savefig(fig, [GroupAnalysisDir saveName '.fig'],'compact');
        saveas(fig, [GroupAnalysisDir saveName '.png']);     
    end
end
close all % close all figure
