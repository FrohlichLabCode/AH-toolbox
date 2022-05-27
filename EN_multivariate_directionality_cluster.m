cluster = 0;
tic

if cluster == 0
    addpath(genpath( 'C:\Users\FrohlichLab\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\'));
    ephysDir      = [ 'C:\Users\FrohlichLab\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\toAnalyze\0153_AttentionTask2Ephys_01_20171020_171020_141018\'];
    saveRootPath = [ 'C:\Users\FrohlichLab\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\Analyzed\0153_AttentionTask2Ephys_01_20171020_171020_141018\conditionalGC\'];
    mkdir(saveRootPath)
    
elseif cluster == 1
    % % Required for use of KilLDevil and parfor
    % ClusterInfo.setQueueName('week')
    % ClusterInfo.setUserDefinedOptions('-M72')
    % ClusterInfo.setUserDefinedOptions('-R mem128')
    % % add paths for cluster computation
    addpath(genpath('/nas/longleaf/home/negahban/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
    pathDir      = [ '/pine/scr/n/e/negahban/Ferret_Grant/toAnalyze/0153_AttentionTask2Ephys_01_20171020_171020_141018']; % CHANGE FOR KILLDEVIL VS LONGLEAF
    saveRootPath = [ '/pine/scr/n/e/negahban/Ferret_Grant/analyzed/0153_AttentionTask2Ephys_01_20171020_171020_141018/conditional_GC/high_freq'];
    mkdir(saveRootPath)
    %code for initialising parallel computing

    numCore = 23; % USR DEFINE
    parpool('local',numCore);
    fileInfo     = dir([pathDir '0*']);
end



%%

%% load data and functions

if cluster == 0
    load('Z:\Ferret Data\0153\BehavioralTraining\VA_training\0153_Level2b_01_20171020.mat')
    cd(ephysDir); % change directory to save data
    addpath('C:\Users\FrohlichLab\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys'); % all the functions needed
    addpath('C:\Users\FrohlichLab\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\mvgc_v1.0')
    % load lfpData
    load([ephysDir 'lfp\lfpMat.mat'])
    % load trigger data
    load([ephysDir 'triggerData.mat'])
    
    
elseif cluster == 1
    % load behavirol data
    load('/pine/scr/n/e/negahban/Ferret_Grant/Behavioral_Data/VA_Training/0153_Level3_11_20171020.mat')
    %load('Z:\Ferret Data\0151\BehavioralTraining\VA_Training\0151_Level6b_01_20170920_behav.mat')
    
    %ephysDir = 'C:\Users\FrohlichLab\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0151\toAnalyze\0151_AttentionTask6b2_1_20170920_170920_181559\';
    ephysDir = '/pine/scr/n/e/negahban/Ferret_Grant/toAnalyze/0153_AttentionTask2Ephys_01_20171020_171020_141018/';
    cd(ephysDir); % change directory to save data
    addpath('/nas/longleaf/home/negahban/Code/'); % all the functions needed
    addpath('/nas/longleaf/home/negahban/Code/mvgc_v1.0/')
    % load lfpData
    load([ephysDir '/lfp/lfpMat.mat'])
    % load trigger data
    load([ephysDir '/triggerData.mat'])
end





%%
% Define frequencies of interest. Linear spacing for Phase slope index, and
% logarithmic spacing for all other methods.
numFreqs = 100;
lowFreq  = 2;
highFreq = 30;
foi      = linspace(lowFreq,highFreq,numFreqs); % linear spacing
% foi   = logspace(log10(lowFreq),log10(highFreq),numFreqs); % log spacing

fs = lfpFs;
twin = [-4,2]; %<<<-------------------------------- median reaction time
onInitInd = 1;
stimTouchInd = 2;

rawFs = 30000;
trialOnset = find(diff(triggerData(onInitInd,:))==1)./rawFs;
trialInit = find(diff(triggerData(onInitInd,:))==-1)./rawFs;
stimOnset = find(diff(triggerData(stimTouchInd,:))==1)./rawFs;
touch = find(diff(triggerData(stimTouchInd,:))==-1)./rawFs;


CalculateMedian = 1;
if CalculateMedian == 1
    % Lowpass filter data at 100Hz
    [b,a] = butter(2,100/(lfpFs/2),'low');
    fdat  = filtfilt(b,a,lfpMat')';
    fdat = lfpMat;
    % For 0153 level 2b
    ppcChans = 1:32;
    pulChans = 41:56;
    vcChans = 65:80;

    % Resample to 200Hz
    idat = median(fdat(ppcChans,:));
    jdat = median(fdat(pulChans,:));
    zdat = median(fdat(vcChans,:));

    lfpMat = [idat;jdat;zdat];
end

region1Chn = 1; %1:32;%9:24;  %<<<------------ PPC
region2Chn = 2; %41:56; %40:56; %33:64; 33:48  %<<<------------ LP/Pulvinar
region3Chn = 3; %65:80; % VC
numChnReg1 = numel(region1Chn);
numChnReg2 = numel(region2Chn);
numChnReg3 = numel(region3Chn);

%lfpMat = lfpPCA;
region1LFP = lfpMat(region1Chn,:);
region2LFP = lfpMat(region2Chn,:);
region3LFP = lfpMat(region3Chn,:);

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

%% within region
% temp6 = 0;
% tvecGC = -1; % initialization of this variable is required before using in parfore
% randReg1 = randperm(32, 10); % PPC
% randReg2 = randperm(16, 8)+40; % LP/Pulvinar
% randReg3 = randperm(16, 8)+64; % VC
for iReg1Chn = 1 %:10%numChnReg1 %or parfor
    xser = region1LFP(iReg1Chn,:);
    temp2 = region2LFP; % this is required for parfor, otherwise it will result a warning saying that overhead communication between workers may happen since the variabel is a broadcast variable.
    temp3 = region3LFP; % this is required for parfor, otherwise it will result a warning saying that overhead communication between workers may happen since the variabel is a broadcast variable.
    temp4 = region2Chn;
    temp5 = region3Chn;
    for iReg2Chn = 1 %:8%numChnReg2
        
        yser = temp2(iReg2Chn,:);
        
        for iReg3Chn = 1 %:8%numChnReg3
            
            zser = temp3(iReg3Chn,:);
            
                    fprintf('\nWorking on pair %d/%d ================================ \n',...
                        iReg3Chn+(iReg2Chn-1)*iReg3Chn+(iReg1Chn-1)*(iReg3Chn+(iReg2Chn-1)*iReg3Chn) , numChnReg1*numChnReg2*numChnReg3);
            
            
            %             % clear spectral data from memory and compute event-triggered power spectrograms
            %             clear xser yser zser
            
            
            
            funcCon = functionalConnectivity_multivariate_cluster(xser,...
                yser,zser,fs,event,twin,region1Chn(iReg1Chn),temp4(iReg2Chn), temp5(iReg3Chn), saveRootPath);
            close all
            
            GC_XtoY_All(iReg1Chn,iReg2Chn,iReg3Chn,:,:) = funcCon.grangerCausality.X_to_Y;
            GC_YtoX_All(iReg1Chn,iReg2Chn,iReg3Chn,:,:) = funcCon.grangerCausality.Y_to_X;
            GC_XtoZ_All(iReg1Chn,iReg2Chn,iReg3Chn,:,:) = funcCon.grangerCausality.X_to_Z;
            GC_ZtoX_All(iReg1Chn,iReg2Chn,iReg3Chn,:,:) = funcCon.grangerCausality.Z_to_X;
            GC_YtoZ_All(iReg1Chn,iReg2Chn,iReg3Chn,:,:) = funcCon.grangerCausality.Y_to_Z;
            GC_ZtoY_All(iReg1Chn,iReg2Chn,iReg3Chn,:,:) = funcCon.grangerCausality.Z_to_Y;
            
            % for inverted time series control
            GC_XtoY_All_invert(iReg1Chn,iReg2Chn,iReg3Chn,:,:) = funcCon.grangerCausality.X_to_Y_invert;
            GC_YtoX_All_invert(iReg1Chn,iReg2Chn,iReg3Chn,:,:) = funcCon.grangerCausality.Y_to_X_invert;
            GC_XtoZ_All_invert(iReg1Chn,iReg2Chn,iReg3Chn,:,:) = funcCon.grangerCausality.X_to_Z_invert;
            GC_ZtoX_All_invert(iReg1Chn,iReg2Chn,iReg3Chn,:,:) = funcCon.grangerCausality.Z_to_X_invert;
            GC_YtoZ_All_invert(iReg1Chn,iReg2Chn,iReg3Chn,:,:) = funcCon.grangerCausality.Y_to_Z_invert;
            GC_ZtoY_All_invert(iReg1Chn,iReg2Chn,iReg3Chn,:,:) = funcCon.grangerCausality.Z_to_Y_invert;
            
%             if iReg1Chn == 1 % This is how I handelled tvecGC in parfor loop
%                 tvecGC(iReg1Chn,iReg2Chn,iReg3Chn,:) = squeeze(funcCon.grangerCausality.tvec);
%             end
            
            
            
        end
    end
    
end
tvecGC =  squeeze(funcCon.grangerCausality.tvec);

%% Eliminate nan sessions
% The results from some sessions have nan elements. Those sessions need to
% be eliminated. Note that simple use of "nanmean" does not work at this stage, but will be used later (EN).
delete_nan = 0;
if delete_nan == 1
    for i=1:size(GC_XtoY_All, 1)
        for j=1:size(GC_XtoY_All, 2)
            for j2 = 1:size(GC_XtoY_All, 3)
                temp1 = squeeze(GC_XtoY_All(i, j, j2, :, :));
                temp2 = isnan(temp1);
                if any(temp2(:))
                    fprintf('\n[i, j]= [%d, %d], NaN in xSpec   ', i, j);
                    GC_XtoY_All(i, j, j2, :, :) = nan;
                end
                clear temp1 temp2
                
                
                temp1 = squeeze(GC_YtoX_All(i, j, j2, :, :));
                temp2 = isnan(temp1);
                if any(temp2(:))
                    fprintf('\n[i, j]= [%d, %d], NaN in xSpec   ', i, j);
                    GC_YtoX_All(i, j, j2, :, :) = nan;
                end
                clear temp1 temp2
                
                
                temp1 = squeeze(GC_XtoZ_All(i, j, j2, :, :));
                temp2 = isnan(temp1);
                if any(temp2(:))
                    fprintf('\n[i, j]= [%d, %d], NaN in xSpec   ', i, j);
                    GC_XtoZ_All(i, j, j2, :, :) = nan;
                end
                clear temp1 temp2
                
                
                temp1 = squeeze(GC_ZtoX_All(i, j, j2, :, :));
                temp2 = isnan(temp1);
                if any(temp2(:))
                    fprintf('\n[i, j]= [%d, %d], NaN in xSpec   ', i, j);
                    GC_ZtoX_All(i, j, j2, :, :) = nan;
                end
                clear temp1 temp2
                
                temp1 = squeeze(GC_ZtoY_All(i, j, j2, :, :));
                temp2 = isnan(temp1);
                if any(temp2(:))
                    fprintf('\n[i, j]= [%d, %d], NaN in xSpec   ', i, j);
                    GC_ZtoY_All(i, j, j2, :, :) = nan;
                end
                clear temp1 temp2
                
                temp1 = squeeze(GC_YtoZ_All(i, j, j2, :, :));
                temp2 = isnan(temp1);
                if any(temp2(:))
                    fprintf('\n[i, j]= [%d, %d], NaN in xSpec   ', i, j);
                    GC_YtoZ_All(i, j, j2, :, :) = nan;
                end
                clear temp1 temp2
                
            end
        end
    end
end
%% Compute the mean values across channel pairs
% tvec = funcCon.tvec;
% tvecGC = funcCon.grangerCausality.tvec;
% psiFreq = funcCon.psiFreq;
%
% avgXSpec = squeeze(mean(mean( xSpecAll ,1),2));
% avgYSpec = squeeze(mean(mean( ySpecAll ,1),2));
% avgPLV = squeeze(mean(mean( plvAll ,1),2));
% avgCoherency = squeeze(mean(mean( coherencyAll ,1),2));
% avgpsiNorm = squeeze(mean(mean( psiNormAll ,1),2));
avgGC_XtoY = squeeze(nanmean(nanmean( nanmean(GC_XtoY_All ,1),2), 3));
avgGC_YtoX = squeeze(nanmean(nanmean( nanmean(GC_YtoX_All ,1),2), 3));
avgGC_XtoZ = squeeze(nanmean(nanmean( nanmean(GC_XtoZ_All ,1),2), 3));
avgGC_ZtoX = squeeze(nanmean(nanmean( nanmean(GC_ZtoX_All ,1),2), 3));
avgGC_YtoZ = squeeze(nanmean(nanmean( nanmean(GC_YtoZ_All ,1),2), 3));
avgGC_ZtoY = squeeze(nanmean(nanmean( nanmean(GC_ZtoY_All ,1),2), 3));

avgGC_XtoY_invert = squeeze(nanmean(nanmean( nanmean(GC_XtoY_All_invert ,1),2), 3));
avgGC_YtoX_invert = squeeze(nanmean(nanmean( nanmean(GC_YtoX_All_invert ,1),2), 3));
avgGC_XtoZ_invert = squeeze(nanmean(nanmean( nanmean(GC_XtoZ_All_invert ,1),2), 3));
avgGC_ZtoX_invert = squeeze(nanmean(nanmean( nanmean(GC_ZtoX_All_invert ,1),2), 3));
avgGC_YtoZ_invert = squeeze(nanmean(nanmean( nanmean(GC_YtoZ_All_invert ,1),2), 3));
avgGC_ZtoY_invert = squeeze(nanmean(nanmean( nanmean(GC_ZtoY_All_invert ,1),2), 3));


cd(saveRootPath)
save('ConditionalGC_avg.mat','avgGC_XtoY','avgGC_YtoX','avgGC_XtoZ','avgGC_ZtoX',...
    'avgGC_YtoZ','avgGC_ZtoY','foi', '-v7.3');
% save('funcCon_avg.mat','avgXSpec','avgYSpec','avgPLV','avgCoherency','avgpsiNorm','avgGC_XtoY','avgGC_YtoX', '-v7.3');
% save('specAll.mat','tvec','foi','xSpecAll','ySpecAll','-v7.3');
% save('plvAll.mat','plvAll','-v7.3');
% save('coherencyAll.mat','coherencyAll','-v7.3');
save('GCAll.mat','tvecGC','GC_XtoY_All','GC_YtoX_All','GC_XtoZ_All','GC_ZtoX_All','GC_YtoZ_All','GC_ZtoY_All','-v7.3');
% save('psiAll.mat','psiAll','psiFreq','psiNormAll','-v7.3');
fprintf('\nDone saving ============================================\n')


%% plotting
for i=1:1
    
    %         avgXSpec = squeeze(nanmean(nanmean( xSpecAll ,1),2));
    %         avgYSpec = squeeze(nanmean(nanmean( ySpecAll ,1),2));
    %         avgPLV = squeeze(nanmean(nanmean( plvAll ,1),2));
    %         avgCoherencey = squeeze(nanmean(nanmean( coherencyAll ,1),2));
    %         avgGC_XtoY = squeeze(nanmean(nanmean( GC_XtoY_All ,1),2));
    %         avgGC_YtoX = squeeze(nanmean(nanmean( GC_YtoX_All ,1),2));
    
    
    % Compute ticks for plotting
% load('GCtvec.mat');
    fig1 = figure('Position',[-1908          38        1438         930]); 
    myclim = 0.4;
    myylim = 30;
    subplot(231);   imagesc(tvecGC,foi,real(avgGC_XtoZ));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: X to Y')
    set(gca,'YDir','normal','TickDir','out')
    cl = colorbar('northoutside'); ylabel(cl,'GC: PPC to VC','FontSize',12)
    caxis([0 myclim]);    %ylim([2 myylim]); % 90+ Hz has saturated values
    
    subplot(232);   imagesc(tvecGC,foi,real(avgGC_ZtoX));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: X to Y')
    set(gca,'YDir','normal','TickDir','out')
    cl = colorbar('northoutside'); ylabel(cl,'GC:  VC to PPC','FontSize',12)
    caxis([0 myclim]);    %ylim([2 myylim]); % 90+ Hz has saturated values
    
    subplot(233);   imagesc(tvecGC,foi,real(avgGC_XtoY));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: X to Y')
    set(gca,'YDir','normal','TickDir','out')
    cl = colorbar('northoutside'); ylabel(cl,'GC: PPC to LPl','FontSize',12)
    caxis([0 myclim]);   % ylim([2 myylim]); % 90+ Hz has saturated values
    
    subplot(234);   imagesc(tvecGC,foi,real(avgGC_YtoX));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: X to Y')
    set(gca,'YDir','normal','TickDir','out')
    cl = colorbar('northoutside'); ylabel(cl,'GC: LPl to PPC','FontSize',12)
    caxis([0 myclim]);   % ylim([2 myylim]); % 90+ Hz has saturated values
    
    subplot(235);   imagesc(tvecGC,foi,real(avgGC_YtoZ));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: X to Y')
    set(gca,'YDir','normal','TickDir','out')
    cl = colorbar('northoutside'); ylabel(cl,'GC: LPl to VC','FontSize',12)
    caxis([0 myclim]);   % ylim([2 myylim]); % 90+ Hz has saturated values
    
    subplot(236);   imagesc(tvecGC,foi,real(avgGC_ZtoY));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: X to Y')
    set(gca,'YDir','normal','TickDir','out')
    cl = colorbar('northoutside'); ylabel(cl,'GC: VC to LPl','FontSize',12)
    caxis([0 myclim]);   % ylim([2 myylim]); % 90+ Hz has saturated values
    
    colormap(jet)
    %plot invert Granger
    fig2 = figure('Position',[-1908          38        1438         930]); 
    
    subplot(231);   imagesc(tvecGC,foi,real(avgGC_XtoZ_invert));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: X to Y')
    set(gca,'YDir','normal','TickDir','out')
    cl = colorbar('northoutside'); ylabel(cl,'inverted GC: PPC to VC','FontSize',12)
    caxis([0 myclim]);    %ylim([2 myylim]); % 90+ Hz has saturated values
    
    subplot(232);   imagesc(tvecGC,foi,real(avgGC_ZtoX_invert));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: X to Y')
    set(gca,'YDir','normal','TickDir','out')
    cl = colorbar('northoutside'); ylabel(cl,'GC:  VC to PPC','FontSize',12)
    caxis([0 myclim]);    %ylim([2 myylim]); % 90+ Hz has saturated values
    
    subplot(233);   imagesc(tvecGC,foi,real(avgGC_XtoY_invert));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: X to Y')
    set(gca,'YDir','normal','TickDir','out')
    cl = colorbar('northoutside'); ylabel(cl,'GC: PPC to LPl','FontSize',12)
    caxis([0 myclim]);   % ylim([2 myylim]); % 90+ Hz has saturated values
    
    subplot(234);   imagesc(tvecGC,foi,real(avgGC_YtoX_invert));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: X to Y')
    set(gca,'YDir','normal','TickDir','out')
    cl = colorbar('northoutside'); ylabel(cl,'GC: LPl to PPC','FontSize',12)
    caxis([0 myclim]);   % ylim([2 myylim]); % 90+ Hz has saturated values
    
    subplot(235);   imagesc(tvecGC,foi,real(avgGC_YtoZ_invert));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: X to Y')
    set(gca,'YDir','normal','TickDir','out')
    cl = colorbar('northoutside'); ylabel(cl,'GC: LPl to VC','FontSize',12)
    caxis([0 myclim]);   % ylim([2 myylim]); % 90+ Hz has saturated values
    
    subplot(236);   imagesc(tvecGC,foi,real(avgGC_ZtoY_invert));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: X to Y')
    set(gca,'YDir','normal','TickDir','out')
    cl = colorbar('northoutside'); ylabel(cl,'GC: VC to LPl','FontSize',12)
    caxis([0 myclim]);   % ylim([2 myylim]); % 90+ Hz has saturated values
    colormap(jet)
end

%%

savefig(fig1, 'avgFuncCon.fig','compact');
saveas(fig1, 'avgFuncCon.png');
savefig(fig2, 'avgFuncCon_invert.fig','compact');
saveas(fig2, 'avgFuncCon_invert.png');

x = toc;
fprintf('time required =%f sec\n', x);

