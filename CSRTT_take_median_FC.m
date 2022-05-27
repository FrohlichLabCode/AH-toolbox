% calculate median from all channels
clear all
skipRec = 0;
doPlot  = 1;
linORlog = 2;


addpath('E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/');
animals = {'0147'}; %{'0168','0169'};
analysisType = 'FC_validChns_150f';

for iAnimal = 1:numel(animals)
    animalCode = animals{iAnimal};
    %PreprocessDir = ['D:/FerretData/' animalCode '/Preprocessed/'];
    AnalysisDir   = ['E:/FerretData/' animalCode '/Analyzed/'];
    %BehavDatDir   = ['D:/FerretData/' animalCode '/behav/'];

    fileInfo = dir([AnalysisDir animalCode '*']); % detect files to load/convert

% loop through each recording
for irec = 1:numel(fileInfo)
    recName = fileInfo(irec).name; %'0168_Opto_010_20180713'
%     splitName   = strsplit(recName,'_');
%     if datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180712', 'InputFormat', 'yyyyMMdd'); continue;end

    recPath = [AnalysisDir recName '/'];
    fileInfo_pair = dir([AnalysisDir recName '/' analysisType '/*']);
    fileInfo_pair(1:2,:)=[]; % first two entries are empty directory
    for ipair = 1:numel(fileInfo_pair)
        rootAnalysisDir = [fileInfo_pair(ipair).folder '/' fileInfo_pair(ipair).name '/'];
        splitName1 = strsplit(fileInfo_pair(ipair).name,'-');
        regionXname = splitName1{1};
        regionYname = splitName1{2};
        
        if exist([rootAnalysisDir 'funcCon_median_Init.mat'],'file') % check if already preprocessed records
            display([recName ' already processed']);
            if skipRec == 1; continue; end % skip record
        end
        display(['processing lfp for rec ' recName]);
        fileInfo_condSpec = dir([rootAnalysisDir '\specAll*' '.mat']);
        fileInfo_condGC = dir([rootAnalysisDir '\GCAll*' '.mat']);
        fileInfo_condPLV = dir([rootAnalysisDir '\PLVAll*' '.mat']);
%         fileInfo_condCoherency = dir([rootAnalysisDir '\coherencyAll*' '.mat']);
%         fileInfo_condPSI = dir([rootAnalysisDir '\psiAll*' '.mat']);
        
        % go through each condition
        for iCond = 1:numel(fileInfo_condSpec)
            splitName = strsplit(fileInfo_condSpec(iCond).name,{'_','.'});
            condName = splitName{2};
            load([rootAnalysisDir fileInfo_condSpec(iCond).name])
            load([rootAnalysisDir fileInfo_condGC(iCond).name]);
            load([rootAnalysisDir fileInfo_condPLV(iCond).name]);
%             load([rootAnalysisDir fileInfo_condCoherency(iCond).name]);
%             load([rootAnalysisDir fileInfo_condPSI(iCond).name]);
            
            
            % calculate median
            avgXSpec     = squeeze(nanmedian(nanmedian( xSpecAll ,1),2)); % change mean to nanmean
            avgYSpec     = squeeze(nanmedian(nanmedian( ySpecAll ,1),2));
            avgPLV       = squeeze(nanmedian(nanmedian( plvAll ,1),2));
            % the following are not priority and takes too long to run
%             avgXNormed   = squeeze(nanmedian(nanmedian( xSpecNormedAll ,1),2));
%             avgYNormed   = squeeze(nanmedian(nanmedian( ySpecNormedAll ,1),2));
%                 tempReal = nanmedian(nanmedian(real(coherencyAll),1),2);
%                 tempImag = 1i*nanmedian(nanmedian(imag(coherencyAll),1),2);
%             avgCoherency = squeeze(tempReal + tempImag); 
%                 tempReal = nanmedian(nanmedian(real(avgImagZ),1),2);
%                 tempImag = 1i*nanmedian(nanmedian(imag(avgImagZ),1),2); 
%             avgImagZ     = squeeze(tempReal + tempImag);         
%             avgpsiNorm   = squeeze(nanmedian(nanmedian( psiNormAll ,1),2));
%             try
%                 tempReal = nanmedian(nanmedian(real(GC_XtoY_All),1),2);
%                 tempImag = 1i*nanmedian(nanmedian(imag(GC_XtoY_All),1),2);     
%             avgGC_XtoY   = squeeze(tempReal + tempImag); 
%                 tempReal = nanmedian(nanmedian(real(GC_YtoX_All),1),2);
%                 tempImag = 1i*nanmedian(nanmedian(imag(GC_YtoX_All),1),2);    
%             avgGC_YtoX   = squeeze(tempReal + tempImag); 
%             catch
%             end

            save([rootAnalysisDir 'funcCon_median_' condName '.mat'],'avgXSpec','avgYSpec','avgPLV', '-v7.3');
%            save(['funcCon_median_' condNames{condID(iCond)} '.mat'],'avgXSpec','avgYSpec','avgXNormed', 'avgYNormed','avgPLV','avgCoherency','avgImagZ','avgpsiNorm', '-v7.3');
            try
            save([rootAnalysisDir 'GC_medain_' condName '.mat'],'tvecGC','avgGC_XtoY','avgGC_YtoX','-v7.3');
            catch
            end
            fprintf('\nDone saving median ============================================\n')
            
            % plot median
                %% plotting
    if doPlot == 1
        %% plot based on median        
        lowFreq = 2; highFreq = 128; numFreqs = 150;
%        psiFreq = foi;
        % Compute ticks for plotting
        if linORlog == 1            
            foi      = linspace(lowFreq,highFreq,numFreqs);
            fois = [2, 5:5:highFreq];
            tickLabel = string(fois); % generate a string array matches fois {"5","10"...}            
%            psitickLabel = string([round(psiFreq(1)) fois(2:end-1) round(psiFreq(end))]); % generate a string array for phase slope index
        elseif linORlog == 2
            foi   = logspace(log10(lowFreq),log10(highFreq),numFreqs);
            fois = 2.^(log2(lowFreq):1:log2(highFreq)); %[2 4 8 12 16 32 64 128];
            tickLabel = string(fois);
%            psitickLabel = string([round(psiFreq(1)) fois(2:end-1) round(psiFreq(end))]);
        end
        for fi = 1:numel(fois)
            [bi,bb] = sort(abs(foi-fois(fi)));
            tickLoc(fi) = bb(1);
%            [bi,bb] = sort(abs(psiFreq-fois(fi)));
%            psitickLoc(fi) = bb(1);
        end

        % Plot
        screensize = get( groot, 'Screensize' );
        fig = figure('Position',[10 50 screensize(3)-150 screensize(4)-150]); %(x,y,width,height) screensize(3)-100

        % plot power spectrum for signal x
        subplot(3,4,1)
        imagesc(tvec,1:numel(foi),pow2db(avgXSpec));
    %    imagesc(tvec,foi,pow2db(avgXSpec));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]'); % title('Signal X power')
    %    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([22 55]);
        cl = colorbar('northoutside'); ylabel(cl,['Power [dB]: ' regionXname],'FontSize',12)

        % plot power spectrum for signal y
        subplot(3,4,2)
        imagesc(tvec,1:numel(foi),pow2db(avgYSpec));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Signal Y power')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([22 55]);
        cl = colorbar('northoutside'); ylabel(cl,['Power [dB]: ' regionYname],'FontSize',12)
        
        try
        subplot(3,4,3)
        imagesc(tvec,1:numel(foi),avgXNormed);
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]'); % title('Signal X power')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel);
        ylim([tickLoc(1) tickLoc(end)]);
        caxis([0 3]);
        cl = colorbar('northoutside'); ylabel(cl,['% power change: ' regionXname],'FontSize',12);

        % plot power spectrum for signal y
        subplot(3,4,4)
        imagesc(tvec,1:numel(foi),avgYNormed);
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Signal Y power')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        ylim([tickLoc(1) tickLoc(end)]);
        caxis([0 3]);
        cl = colorbar('northoutside'); ylabel(cl,['% power change:' regionYname],'FontSize',12);
        catch
        end
        
        try
        % plot phase locking value
        subplot(3,4,5)
        imagesc(tvec,1:numel(foi),avgPLV);
        %imagesc(tvec,foi,avgPLV);
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([0.1 0.7]);
        cl = colorbar('northoutside'); ylabel(cl,'Phase locking value','FontSize',12)
        catch
        end

        try
        % plot coherence
        subplot(3,4,6)
        imagesc(tvec,1:numel(foi),abs(avgCoherency));
        %imagesc(tvec,foi,abs(avgCoherencey));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Coherence')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        %caxis([0 0.6]);
        cl = colorbar('northoutside'); ylabel(cl,'Coherence','FontSize',12)
        catch
        end

        try
        % plot imaginary coherence
        subplot(3,4,7)
        imagesc(tvec,1:numel(foi),abs(avgImagZ));
        %imagesc(tvec,foi,imag(avgCoherency));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Imaginary coherence')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([0 4]);
        cl = colorbar('northoutside'); ylabel(cl,'Imag coherence (z)','FontSize',12)
        % plot phase slope index
        catch
        end

        try
            subplot(3,4,8)
            imagesc(tvec,1:numel(psiFreq),avgpsiNorm);
            %imagesc(tvec,psiFreq,avgpsiNorm);
            xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Phase slope index')
            ylim([psitickLoc(1) psitickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',psitickLoc,'YTickLabel',psitickLabel)
            %cl = colorbar('northoutside'); ylabel(cl,'Phase slope index (z-score)','FontSize',15)
            cl = colorbar('northoutside'); ylabel(cl,'Phase slope index (z)','FontSize',12)
            % plot granger causality X to Y
            caxis([-4 4])
        catch
        end
        try
            subplot(3,4,9)
            imagesc(tvecGC,1:numel(foi),real(avgGC_XtoY));
            %imagesc(tvecGC,foi,real(avgGC_XtoY));
            xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('GC: X to Y')
            set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
            %cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: X to Y','FontSize',15)
            cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionXname ' to ' regionYname],'FontSize',12)
            caxis([0 0.3]); ylim([tickLoc(1) tickLoc(end)]); % 90+ Hz has saturated values
            % plot granger causality Y to X
            subplot(3,4,10)
            imagesc(tvecGC,1:numel(foi),real(avgGC_YtoX));
            %imagesc(tvecGC,foi,real(avgGC_YtoX));
            xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('GC: Y to X')
            set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
            %cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: Y to X','FontSize',15)
            cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionYname ' to ' regionXname],'FontSize',12)
            caxis([0 0.3]); ylim([tickLoc(1) tickLoc(end)]) % 90+ Hz has saturated values
        catch
        end
        colormap(jet)

        savefig(fig, [rootAnalysisDir 'medianFuncCon_' condName '_' num2str(lowFreq) '-' num2str(highFreq) 'Hz.fig'],'compact');
        saveas(fig, [rootAnalysisDir 'medianFuncCon_' condName '_' num2str(lowFreq) '-' num2str(highFreq) 'Hz.png']);
        close all
        clear avgXSpec avgYSpec avgXNormed avgYNormed
        clear avgPLV avgCoherency avgImagZ avgpsiNorm avgGC_XtoY avgGC_YtoX
    end
        end
    end
end
end

            