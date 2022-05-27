clear

cluster = 0;
tic

%% load data and functions

if cluster == 0
    addpath('J:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys'); % all the functions needed
    addpath('J:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\mvgc_v1.0') % for Granger Causality
    ephysDir     = ['J:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\toAnalyze\0153_AttentionTask3_13\'];     
    saveRootPath = ['J:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\Analyzed\0153_AttentionTask3_13\1-80Hzlog\']; 
    mkdir(saveRootPath);  
    %addpath(genpath(ephysDir)); % add folder and subfolders to path
    cd(ephysDir); % change directory to save data

    % load behavioral data
    load('J:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\5CSRTT\Behav data\0153\0153_Level3_13_20171025.mat'); %20,24,25,27
    % load lfpData
    load([ephysDir 'lfp\lfpMat.mat']);
    % load trigger data
    load([ephysDir 'triggerData.mat']);
        
elseif cluster == 1
    % % Required for use of KilLDevil and parfor
    % ClusterInfo.setQueueName('week')
    % ClusterInfo.setUserDefinedOptions('-M72')
    % ClusterInfo.setUserDefinedOptions('-R mem128')
    % % add paths for cluster computation
    addpath(genpath('/nas/longleaf/home/angelvv/Code/')); % all the functions needed
    % CHANGE FOR KILLDEVIL VS LONGLEAF
    ephysDir     = [ '/pine/scr/h/w/angelvv/5CSRTT/toAnalyze/0147_AttentionTask6_test_20170824']; % CHANGE FOR KILLDEVIL VS LONGLEAF
    saveRootPath = [ '/pine/scr/h/w/angelvv/5CSRTT/Analyzed/0147_AttentionTask6_test_20170824/LPl_PPC_touch'];
    mkdir(saveRootPath)
    cd(ephysDir)
    
    % load behavirol data
    load('/pine/scr/h/w/angelvv/5CSRTT/Behavioral/0147_Level6_01_20170824_1_behav.mat')
    % load lfpData
    load([ephysDir '/lfp/lfpMat.mat'])
    % load trigger data
    load([ephysDir '/triggerData.mat'])

    %code for initialising parallel computing    
    numCore = 24; % USER DEFINE
    parpool('local',numCore);
    fileInfo = dir([pathDir '0*']);
end

%%
% Define frequencies of interest. Linear spacing for Phase slope index, and spacing for all other methods.
linORlog = 2;
if linORlog == 1
    numFreqs = 200;
    lowFreq  = 1;
    highFreq = 80;
    foi      = linspace(lowFreq,highFreq,numFreqs); % linear spacing
end
if linORlog == 2
    numFreqs = 200;
    lowFreq  = 1;
    highFreq = 80;
    foi   = logspace(log10(lowFreq),log10(highFreq),numFreqs); % log spacing
end


fs = lfpFs;
twin = [-4 2]; %<<<-------------------------------- median reaction time
onInitInd = 1;
stimTouchInd = 2;

rawFs = 30000;
trialOnset = find(diff(triggerData(onInitInd,:))==1)./rawFs;
trialInit = find(diff(triggerData(onInitInd,:))==-1)./rawFs;
stimOnset = find(diff(triggerData(stimTouchInd,:))==1)./rawFs;
touch = find(diff(triggerData(stimTouchInd,:))==-1)./rawFs;

%% label channels
doPCA = 1; channelMap = 2;
if channelMap == 1 
    validChn{1} = 1:32;  %<<<------------ PPC
    validChn{2} = 41:56; %<<<------------ LPl
    validChn{3} = 65:80; %<<<------------ VC
end
if channelMap == 2
    validChn{1} = 1:32; %1:32;%9:24;  %<<<------------ PPC
    validChn{2} = 33:48; %40:56; %33:64; 33:48  %<<<------------ LPl
    validChn{3} = 49:64; % VC
end
if channelMap == 3
    validChn{1} = 1:32; %1:32;%9:24;  %<<<------------ Pulvinar
    validChn{2} = 33:64; %40:56; %33:64; 33:48  %<<<------------ PPC
    validChn{3} = 65:80; % VC
end

% manually exclude noisy channels from 3 regions
deleteChn = [5,9,10,11,12,20,23,38,43,50:56,58]; %[9,10,11,12,20,43,50:56];
tic
for iRegion = 1:3 % 5sec
    validChn{iRegion} = setdiff(validChn{iRegion},deleteChn);
    nValidChn(iRegion) = numel(validChn{iRegion});
    
    if doPCA == 1
        [coeff{iRegion},score{iRegion},latent{iRegion},tsquared{iRegion},explained{iRegion},mu{iRegion}] ...
        	= pca(lfpMat(validChn{iRegion},:)', 'VariableWeights','variance', 'Centered',true); % same result
        lfpMedian(iRegion,:) = median(lfpMat(validChn{iRegion},:),1); % median lfp across channels
        lfpPCA(iRegion,:) = mean(lfpMat(validChn{iRegion},:).*repmat(coeff{iRegion}(:,1),[1 length(lfpMat)]),1); % get PCA1 for each region (weighted avg of each channel's lfp, weight=PCA1 coefficient)           
        
        % scale up lfpPCA
        PCAamp = mean(abs(lfpPCA(iRegion,:))); %VCamp=7
        lfpamp = mean(mean(abs(lfpMat(validChn{iRegion},:)),2),1);%VCamp=29
        lfpPCA(iRegion,:) = lfpPCA(iRegion,:)*lfpamp/PCAamp;%VCamp=29
    end
end

toc


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

event = touch(hitTrials);
event = event(2:length(event)); % exam the trials and note to see whether the animal performs normally

save([saveRootPath 'lfpPCA.mat'],'lfpFs','lfpPCA','validChn','event','twin','-v7.3');

%% within region
doFuncCon = 0;
if doFuncCon == 1
for iReg1Chn = 1:2 
    
    xser = lfpPCA(iReg1Chn,:);
    
    for iReg2Chn = 2:3
%         fprintf('\nWorking on pair %d/%d, OuterCounter=%d, InnerCounter=%d ================================ \n',...
%             iReg2Chn+(iReg1Chn-1)*numChnReg3 , numChnReg1*numChnReg3, iReg1Chn, iReg2Chn);
        if iReg1Chn == iReg2Chn % skip this condition
            continue;
        end
        yser = lfpPCA(iReg2Chn,:);
        funcCon = is_functionalConnectivity(xser,yser,fs,event,twin,iReg1Chn,iReg2Chn, saveRootPath);
        close all
        
        % for each struct matrix, it's frequency by time
        xSpecAll(iReg1Chn,iReg2Chn,:,:) = funcCon.xspec;
        ySpecAll(iReg1Chn,iReg2Chn,:,:) = funcCon.yspec;
        plvAll(iReg1Chn,iReg2Chn,:,:) = funcCon.plv;
        coherencyAll(iReg1Chn,iReg2Chn,:,:) = funcCon.coherency; % might want to keep the complex part to look at phase lags
        coherenceAll(iReg1Chn,iReg2Chn,:,:) = funcCon.coherence;
        iCoherenceAll(iReg1Chn,iReg2Chn,:,:) = funcCon.imaginaryCoherence;
        imagZAll(iReg1Chn,iReg2Chn,:,:) = funcCon.imagZ;
        psiAll(iReg1Chn,iReg2Chn,:,:) = funcCon.psi;
        psiNormAll(iReg1Chn,iReg2Chn,:,:) = funcCon.psiNorm;
        GC_XtoY_All(iReg1Chn,iReg2Chn,:,:) = funcCon.grangerCausality.X_to_Y;
        GC_YtoX_All(iReg1Chn,iReg2Chn,:,:) = funcCon.grangerCausality.Y_to_X;
        
        
        cd(saveRootPath)
        save(['funcConPCA_region' num2str(iReg1Chn) '-' num2str(iReg2Chn) '.mat'],'funcCon','foi','event','twin','-v7.3');
        
    end
    
end
end

%% Eliminate nan sessions
% % The results from some sessions have nan elements. Those sessions need to
% % be eliminated. Note that simple use of "nanmean" does not work at this stage, but will be used later (EN).
% 
% for i=1:size(xSpecAll, 1)
%     for j=1:size(xSpecAll, 2)
%         temp1 = squeeze(xSpecAll(i, j, :, :));
%         temp2 = isnan(temp1);
%         if any(temp2(:))
%             fprintf('\n[i, j]= [%d, %d], NaN in xSpec   ', i, j);
%             xSpecAll(i, j, :, :) = nan;
%         end
%         clear temp1 temp2
%         
%         
%         temp1 = squeeze(ySpecAll(i, j, :, :));
%         temp2 = isnan(temp1);
%         if any(temp2(:))
%             fprintf('\n[i, j]= [%d, %d], NaN in ySpec   ', i, j);
%             ySpecAll(i, j, :, :) = nan;
%         end
%         clear temp1 temp2
%         
%         temp1 = squeeze(GC_XtoY_All(i, j, :, :));
%         temp2 = isnan(temp1);
%         if any(temp2(:))
%             fprintf('\n[i, j]= [%d, %d], NaN in GC_XtoY   ', i, j);
%             GC_XtoY_All(i, j, :, :) = nan;
%         end
%         clear temp1 temp2
%         
%         temp1 = squeeze(GC_YtoX_All(i, j, :, :));
%         temp2 = isnan(temp1);
%         if any(temp2(:))
%             fprintf('\n[i, j]= [%d, %d], NaN in GC_YtoX   ', i, j);
%             GC_YtoX_All(i, j, :, :) = nan;
%         end
%         clear temp1 temp2
%         
%         temp1 = squeeze(plvAll(i, j, :, :));
%         temp2 = isnan(temp1);
%         if any(temp2(:))
%             fprintf('\n[i, j]= [%d, %d], NaN in plv   ', i, j);
%             plvAll(i, j, :, :) = nan;
%         end
%         clear temp1 temp2
%         
%         temp1 = squeeze(coherencyAll(i, j, :, :));
%         temp2 = isnan(temp1);
%         if any(temp2(:))
%             fprintf('\n[i, j]= [%d, %d], NaN in coherency   ', i, j);
%             coherencyAll(i, j, :, :) = nan;
%         end
%         clear temp1 temp2
%         
%         temp1 = squeeze(imagZAll(i, j, :, :));
%         temp2 = isnan(temp1);
%         if any(temp2(:))
%             fprintf('\n[i, j]= [%d, %d], NaN in coherency   ', i, j);
%             imagZAll(i, j, :, :) = nan;
%         end
%         clear temp1 temp2
%         
%         
%     end
% end

%% Compute the mean values across channel pairs
% tvec = funcCon.tvec;
% tvecGC = funcCon.grangerCausality.tvec;
% psiFreq = funcCon.psiFreq;
% 
% cd(saveRootPath)
% %save('funcCon_avg.mat','avgXSpec','avgYSpec','avgPLV','avgCoherency','avgImagZ','avgpsiNorm', '-v7.3');
% save('specAll.mat','tvec','foi','xSpecAll','ySpecAll','-v7.3');
% save('plvAll.mat','plvAll','-v7.3');
% save('coherencyAll.mat','coherencyAll','-v7.3');
% save('GCAll.mat','tvecGC','GC_XtoY_All','GC_YtoX_All','-v7.3');
% save('psiAll.mat','psiAll','psiFreq','psiNormAll','-v7.3');
% fprintf('\nDone saving ============================================\n')


%% plotting
doPlot = 0;
if doPlot == 1
    
    avgXSpec = squeeze(nanmean(nanmean( xSpecAll ,1),2));
    avgYSpec = squeeze(nanmean(nanmean( ySpecAll ,1),2));
    avgPLV = squeeze(nanmean(nanmean( plvAll ,1),2));
    avgCoherency = squeeze(nanmean(nanmean( coherencyAll ,1),2));
    avgGC_XtoY = squeeze(nanmean(nanmean( GC_XtoY_All ,1),2));
    avgGC_YtoX = squeeze(nanmean(nanmean( GC_YtoX_All ,1),2));

    
    % Compute ticks for plotting
    screensize = get( groot, 'Screensize' );
    fig = figure('Position',[10 50 screensize(3)-100 screensize(4)-150]);
    
    % Compute ticks for plotting
    % fois = [0.5 1 2 4 8 16 32 64 128];
    fois = [5 10 15 20 25 30];
    for fi = 1:numel(fois)
        [bi,bb] = sort(abs(foi-fois(fi)));
        tickLoc(fi) = bb(1);
    end
    tickLabel = string(fois);
    
    % do the same for phase slope index frequencies
    % fois = [0.6 1 2 4 8 16 32 64];
    fois = [5 10 15 20 25 30];
    for fi = 1:numel(fois)
        [bi,bb] = sort(abs(funcCon.psiFreq-fois(fi)));
        psitickLoc(fi) = bb(1);
    end
    psitickLabel = string(fois);
    
    try
    % plot power spectrum for signal x
    subplot(2,4,1)
    imagesc(tvec,1:numel(foi),pow2db(avgXSpec));
%    imagesc(tvec,foi,pow2db(avgXSpec));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)'); % title('Signal X power')
%    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    ylim([2 60]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    caxis([15 60]);
    cl = colorbar('northoutside'); ylabel(cl,'Power (dB) signal X','FontSize',12)
        catch
    end
    try
    % plot power spectrum for signal y
    subplot(2,4,2)
    imagesc(tvec,1:numel(foi),pow2db(avgYSpec));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('Signal Y power')
    ylim([2 60]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    caxis([15 60]);
    cl = colorbar('northoutside'); ylabel(cl,'Power (dB) signal Y','FontSize',12)
    catch
    end
    try
    % plot phase locking value
    subplot(2,4,3)
    imagesc(tvec,1:numel(foi),avgPLV);
    %imagesc(tvec,foi,avgPLV);
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('PLV')
    ylim([2 60]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    %caxis([0 0.8]);
    cl = colorbar('northoutside'); ylabel(cl,'Phase locking value','FontSize',12)
    catch
    end
    
    try
        % plot coherence
    subplot(2,4,4)
    imagesc(tvec,1:numel(foi),abs(avgCoherency));
    %imagesc(tvec,foi,abs(avgCoherencey));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('Coherence')
    ylim([2 60]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    %caxis([0 0.6]);
    cl = colorbar('northoutside'); ylabel(cl,'Coherence','FontSize',12)
   
    catch
    end
    
    try
        % plot imaginary coherence
    subplot(2,4,5)
    imagesc(tvec,1:numel(foi),abs(avgImagZ));
    %imagesc(tvec,foi,imag(avgCoherency));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('Imaginary coherence')
    ylim([2 60]);caxis([0 0.6]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    %cl = colorbar('northoutside'); ylabel(cl,'Imaginary coherence (z-score)','FontSize',15)
    cl = colorbar('northoutside'); ylabel(cl,'Imag coherence (z)','FontSize',12)
    % plot phase slope index
    catch
    end
    
    try
        subplot(2,4,6)
    imagesc(tvec,1:numel(psiFreq),avgpsiNorm);
    %imagesc(tvec,psiFreq,avgpsiNorm);
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('Phase slope index')
    ylim([2 60]);set(gca,'YDir','normal','TickDir','out','YTick',psitickLoc,'YTickLabel',psitickLabel)
    %cl = colorbar('northoutside'); ylabel(cl,'Phase slope index (z-score)','FontSize',15)
    cl = colorbar('northoutside'); ylabel(cl,'Phase slope index (z)','FontSize',12)
    % plot granger causality X to Y
    catch
    end
    try
        subplot(2,4,7)
        imagesc(tvecGC,1:numel(foi),real(avgGC_XtoY));
        %imagesc(tvecGC,foi,real(avgGC_XtoY));
        xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: X to Y')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        %cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: X to Y','FontSize',15)
        cl = colorbar('northoutside'); ylabel(cl,'GC: X to Y','FontSize',12)
        caxis([0 0.3]); ylim([2 60]); % 90+ Hz has saturated values
        % plot granger causality Y to X
        subplot(2,4,8)
        imagesc(tvecGC,1:numel(foi),real(avgGC_YtoX));
        %imagesc(tvecGC,foi,real(avgGC_YtoX));
        xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: Y to X')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        %cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: Y to X','FontSize',15)
        cl = colorbar('northoutside'); ylabel(cl,'GC: Y to X','FontSize',12)
        caxis([0 0.3]); ylim([2 60]) % 90+ Hz has saturated values
    catch
    end

    colormap(jet)
    
    savefig(fig, 'avgFuncCon.fig','compact');
    saveas(fig, 'avgFuncCon.png');
end

x = toc;
fprintf('time required =%f sec\n', x);

