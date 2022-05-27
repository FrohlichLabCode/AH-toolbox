
%% load data and functions

% load behavirol data
load('Z:\Ferret Data\0147\BehavioralTraining\VA_Training\0147_Level6b_01_20170825_1_behav.mat')
%load('Z:\Ferret Data\0151\BehavioralTraining\VA_Training\0151_Level6b_01_20170920_behav.mat')

%ephysDir = 'C:\Users\FrohlichLab\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0151\toAnalyze\0151_AttentionTask6b2_1_20170920_170920_181559\';
ephysDir = 'F:\DB\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0147\toAnalyze\0147_AttentionTask6b_01_20170825_2_170825_160801\';
cd(ephysDir); % change directory to save data
addpath('F:\DB\Dropbox (Frohlich Lab)\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\'); % all the functions needed --Angel
addpath('F:\DB\Dropbox (Frohlich Lab)\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\mvgc_v1.0')
%addpath('C:\Users\FrohlichLab\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\');
%addpath('C:\Users\FrohlichLab\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\mvgc_v1.0');

% load lfpData
load([ephysDir 'lfp\lfpMat.mat'])
% load trigger data
load([ephysDir 'triggerData.mat'])


%%
% Define frequencies of interest. Linear spacing for Phase slope index, and
% logarithmic spacing for all other methods.
numFreqs = 100;
lowFreq  = 0.5;
highFreq = 128;
% foi      = linspace(lowFreq,highFreq,numFreqs); % linear spacing
foi   = logspace(log10(lowFreq),log10(highFreq),numFreqs); % log spacing

fs = lfpFs;
twin = [-6 2]; %<<<--------------------------------
onInitInd = 1;
stimTouchInd = 2;

rawFs = 20000;
trialOnset = find(diff(triggerData(onInitInd,:))==1)./rawFs;
trialInit = find(diff(triggerData(onInitInd,:))==-1)./rawFs;
stimOnset = find(diff(triggerData(stimTouchInd,:))==1)./rawFs;
touch = find(diff(triggerData(stimTouchInd,:))==-1)./rawFs;

region1Chn = 5:32;%9:24;  %<<<------------
region2Chn = 33:48; %40:56; %33:64; 33:48  %<<<------------
numChnReg1 = numel(region1Chn);
numChnReg2 = numel(region2Chn);

region1LFP = lfpMat(region1Chn,:);
region2LFP = lfpMat(region2Chn,:);

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

event = stimOnset(hitTrials);

%% within region 

for iReg1Chn = 1:numChnReg1
    iReg1Chn

    xser = region1LFP(iReg1Chn,:);
    
    for iReg2Chn = 1:numChnReg2
        
        yser = region2LFP(iReg2Chn,:);
        
        
        funcCon = is_functionalConnectivity(xser,yser,fs,event,twin,region1Chn(iReg1Chn),region2Chn(iReg2Chn));
        
        close all
        
        % for each struct matrix, it's frequency by time
        xSpecAll(iReg1Chn,iReg2Chn,:,:) = funcCon.xspec;
        ySpecAll(iReg1Chn,iReg2Chn,:,:) = funcCon.yspec;  
        %only need to change iReg2Chn
        %but Subscripted assignment dimension mismatch.for
        %if iReg1Chn == 1;  ySpecAll(iReg2Chn,:,:) = funcCon.yspec;
        plvAll(iReg1Chn,iReg2Chn,:,:) = funcCon.plv;
        coherencyAll(iReg1Chn,iReg2Chn,:,:) = funcCon.coherency; % might want to keep the complex part to look at phase lags 
        coherenceAll(iReg1Chn,iReg2Chn,:,:) = funcCon.coherence;
        iCoherenceAll(iReg1Chn,iReg2Chn,:,:) = funcCon.imaginaryCoherence;
        imagZAll(iReg1Chn,iReg2Chn,:,:) = funcCon.imagZ;
        psiAll(iReg1Chn,iReg2Chn,:,:) = funcCon.psi;
        psiNormAll(iReg1Chn,iReg2Chn,:,:) = funcCon.psiNorm;
        GC_XtoY_All(iReg1Chn,iReg2Chn,:,:) = funcCon.grangerCausality.X_to_Y;
        GC_YtoX_All(iReg1Chn,iReg2Chn,:,:) = funcCon.grangerCausality.Y_to_X;
        
        
    end
   
end

tvec = funcCon.tvec;
tvecGC = funcCon.grangerCausality.tvec;
psiFreq = funcCon.psiFreq;

avgXSpec = squeeze(mean(mean( xSpecAll ,1),2));
avgYSpec = squeeze(mean(mean( ySpecAll ,1),2));
avgPLV = squeeze(mean(mean( plvAll ,1),2));
avgCoherence = squeeze(mean(mean( coherenceAll ,1),2));
avgimagZ = squeeze(mean(mean(abs( imagZAll) ,1),2));
avgpsiNorm = squeeze(mean(mean( psiNormAll ,1),2));

avgGC_XtoY = squeeze(mean(mean( GC_XtoY_All ,1),2));
avgGC_YtoX = squeeze(mean(mean( GC_YtoX_All ,1),2));

save('funcCon_avg.mat','avgXSpec','avgYSpec','avgPLV','avgCoherence','avgimagZ','avgpsiNorm','avgGC_XtoY','avgGC_YtoX','-v7.3');
save('specAll.mat','tvec','foi','xSpecAll','ySpecAll','-v7.3');
save('plvAll.mat','plvAll','-v7.3');
save('coherencyAll.mat','coherencyAll','-v7.3');
% save('coherenceAll.mat','coherenceAll','-v7.3'); % too big to save
% save('iCoherenceAll','iCoherenceAll','-v7.3');
% save('imagZAll.mat','imagZAll','-v7.3');
% save('psiAll.mat','psiAll','psiFreq','psiNormAll','-v7.3');
% save('GCAll.mat','tvecGC','GC_XtoY_All','GC_YtoX_All','-v7.3');

%% plotting

% Compute ticks for plotting
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 screensize(3)-100 screensize(4)-150]);
% Compute ticks for plotting
fois = [0.5 1 2 4 8 16 32 64 128];
for fi = 1:numel(fois)
    [bi,bb] = sort(abs(foi-fois(fi)));
    tickLoc(fi) = bb(1);
end
tickLabel = {'0.5','1','2','4','8','16','32','64','128'};

% do the same for phase slope index frequencies
fois = [0.6 1 2 4 8 16 32 64];
for fi = 1:numel(fois)
    [bi,bb] = sort(abs(funcCon.psiFreq-fois(fi)));
    psitickLoc(fi) = bb(1);
end
psitickLabel = {'0.6','1','2','4','8','16','32','64'};

% plot power spectrum for signal x
subplot(2,4,1)
imagesc(tvec,1:numel(foi),pow2db(avgXSpec));
xlabel('Time to event (s)'); ylabel('Frequency (Hz)'); % title('Signal X power')
set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
cl = colorbar('northoutside'); ylabel(cl,'Power (dB) signal X','FontSize',12)
% plot power spectrum for signal y
subplot(2,4,2)
imagesc(tvec,1:numel(foi),pow2db(avgYSpec));
xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('Signal Y power')
set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
cl = colorbar('northoutside'); ylabel(cl,'Power (dB) signal Y','FontSize',12)
% plot phase locking value
subplot(2,4,3)
imagesc(tvec,1:numel(foi),avgPLV);
xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('PLV')
set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
cl = colorbar('northoutside'); ylabel(cl,'Phase locking value','FontSize',12)
% plot coherence
subplot(2,4,4)
imagesc(tvec,1:numel(foi),avgCoherence);
xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('Coherence')
set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
cl = colorbar('northoutside'); ylabel(cl,'Coherence','FontSize',12)
% plot imaginary coherence
subplot(2,4,5)
imagesc(tvec,1:numel(foi),avgimagZ);
xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('Imaginary coherence')
set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
%cl = colorbar('northoutside'); ylabel(cl,'Imaginary coherence (z-score)','FontSize',15)
cl = colorbar('northoutside'); ylabel(cl,'Imag coherence (z)','FontSize',12)
% plot phase slope index
subplot(2,4,6)
imagesc(tvec,1:numel(psiFreq),avgpsiNorm);
xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('Phase slope index')
set(gca,'YDir','normal','TickDir','out','YTick',psitickLoc,'YTickLabel',psitickLabel)
%cl = colorbar('northoutside'); ylabel(cl,'Phase slope index (z-score)','FontSize',15)
cl = colorbar('northoutside'); ylabel(cl,'Phase slope index (z)','FontSize',12)
% plot granger causality X to Y
try
    subplot(2,4,7)
    imagesc(tvecGC,1:numel(foi),avgGC_XtoY);
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: X to Y')
    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    %cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: X to Y','FontSize',15)
    cl = colorbar('northoutside'); ylabel(cl,'GC: X to Y','FontSize',12)
    caxis([0 0.3]); ylim([1 90]); % 90+ Hz has saturated values
    % plot granger causality Y to X
    subplot(2,4,8)
    imagesc(tvecGC,1:numel(foi),avgGC_YtoX);
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: Y to X')
    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    %cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: Y to X','FontSize',15)
    cl = colorbar('northoutside'); ylabel(cl,'GC: Y to X','FontSize',12)
    caxis([0 0.3]); ylim([1 90]) % 90+ Hz has saturated values
catch
end
colormap(jet)

savefig(fig, 'avgFuncCon.fig','compact');
saveas(fig, 'avgFuncCon.png');