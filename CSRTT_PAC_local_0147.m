% This code is adopted from Iain's working code, the result doesn't match
% what we see in the PAC_PLV_method, which is the conventional way of
% calculating PAC and also the way we used in our paper.
% So be careful when using this code. 


clear

cluster = 0;
tic

animalCodes = {'0171'};
level = '6b'; % 
skipRec = 1;
MedianorPCA = 3; %0=_validChns, 1=mdChn, 2=PCA, 3=opto1Chn, 4=_validAnaChns
analysisType = 'PAC';
alignID = 2; %1=Init, 2=Stim, 3=Touch, 4=Opto
hitMissID = 1; %1=Correct, 2=Premature, 3=Incorrect, 4=Omission, 5=noPremature

%% Define frequencies of interest.
[foi, tickLoc, tickLabel,psitickLoc,psitickLabel] = getFoiLabel(2, 128, 150, 2); % (lowFreq, highFreq, numFreqs, linORlog)
numFreqs     = numel(foi);
lfpFs        = 1000;

%% load data and functions
for iAnimal = 4%1:numel(animals)
    animalCode = animals{iAnimal};
    if strcmp(animalCode,'0171') && level(1) == '6'
        doMix = 1;
    else
        doMix = 0;
    end
    if doMix == 1; mixSuffix = '_mix'; else; mixSuffix = []; end

if cluster == 0 %linux use '/', windows matlab can use both '/' and '\'
    addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
    addpath(genpath('E:\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\Toolboxes\eeglab2019_0\'));
    PreprocessDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/Preprocessed' mixSuffix '/'];
    AnalysisDir   = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/Analyzed/'];
    GroupAnalysisDir = ['E:/FerretData/' animalCode '/GroupAnalysis/' analysisType folderSuffix '_' level '/'];
elseif cluster == 1
    addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
    PreprocessDir = ['/pine/scr/a/n/angelvv/FerretData/' animalCode '/Preprocessed' mixSuffix '/']; % CHANGE FOR KILLDEVIL VS LONGLEAF
    AnalysisDir   = ['/pine/scr/a/n/angelvv/FerretData/' animalCode '/Analyzed/'];
    GroupAnalysisDir = ['/pine/scr/a/n/angelvv/FerretData/' animalCode '/GroupAnalysis/' analysisType folderSuffix '_' level '/'];
    
    %code for initialising parallel computing
    if doParfor == 1
    numCore = 16; % USR DEFINE, max 24 physical + 24 virtual core per computer
    myPool = parpool('local',numCore,'SpmdEnabled',false);  
    end
end

[alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level);
region = getAnimalInfo(animalCode);
alignName = alignNames{alignID}; %Init
hitMissName = hitMissNames{hitMissID}(1:3); %Cor
alignHitName = [alignName hitMissName]; %InitCorAll

if level(1) == '7' || level(1) == '8'
    [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'optoEventTimes_' alignHitName '.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
    condIDs = [1,2,5]; %[1,2,3,4,5] 1=Theta,2=Alpha,3=ArTheta,4=ArAlpha,5=Sham
elseif level(1) == '9'
    [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'optoEventTimes_' alignHitName '.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
    condIDs = [2,5]; %[1,2,3,4,5] 1=Theta,2=Alpha,3=ArTheta,4=ArAlpha,5=Sham    
elseif level(1) == '6'
    [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'eventTimes_' alignHitName '.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
    condIDs = [4];%1=D4, 2=D5,3=D6,4=Dall    
end

region = getAnimalInfo(animalCode);
folderSuffix = getFolderSuffix(MedianorPCA); %0=_validChns; 1=_median; 2=_PCA; 3=_firstChn;
fileInfo   = dir([PreprocessDir animalCode '_Level' level '*']); % detect fileInfo to load/convert  '_LateralVideo*' can't process opto

%% Load each session
for irec = 1:numel(fileInfo)     
    recName = fileInfo(irec).name;
    splitName   = strsplit(recName,'_');
    sessionID   = splitName{3};
    rootPreprocessDir = [PreprocessDir recName '/'];
    rootAnalysisDir   = [AnalysisDir recName '/' analysisType folderSuffix '/'];

    %% load lfp
    fprintf('\nWorking on record %s =============== \n',recName');
    [regionChn, regionLFP, ~, lfpFs] = getRegionLFP(rootPreprocessDir, MedianorPCA, newFs);

    % get LPl and VC from:
%load('E:\FerretData\0147\Preprocessed\0147_AttentionTask6a_03_20170926\lfp\lfpMat.mat')
load('E:\FerretData\0147\Preprocessed\0147_AttentionTask6a_04_20170927\lfp\lfpMat.mat')
LPl = lfpMat(33:48,2*60*lfpFs:7*60*lfpFs); %17min
VC  = lfpMat(1:32,2*60*lfpFs:7*60*lfpFs);
% get PPC from:
%load('E:\FerretData\0147\Preprocessed\0147_AttentionTask6b_11_20170912\lfp\lfpMat.mat')
load('E:\FerretData\0147\Preprocessed\0147_AttentionTask6b_14_20170915\lfp\lfpMat.mat')
PPC = lfpMat(1:32,2*60*lfpFs:7*60*lfpFs); 


% downsample to 200Hz
newFs   = 200;
[nr,dr] = rat(lfpFs/newFs); %returns two arrays, N and D, such that N./D ~ X

pulLFP = nanmedian(resample(LPl',dr,nr)',1);
ppcLFP = nanmedian(resample(PPC',dr,nr)',1); % not the same length
vcLFP  = nanmedian(resample(VC',dr,nr)',1); 

cFS  = 2; % every cFSth sample to get a sample rate of 100Hz
dvec = 1:cFS:numel(pulLFP);
ppcMat = nan(numFreqs,numFreqs,numel(dvec));
pulMat = nan(numFreqs,numFreqs,numel(dvec));
vcMat  = nan(numFreqs,numFreqs,numel(dvec));
downFS = newFs/cFS;

wavs   = is_makeWavelet(foi,newFs); % each wavs is 1 by time window for 1 cycle at newFs sampling freq

for fi = 1:numFreqs
    display(num2str(fi))
    tmp = conv(ppcLFP,wavs{fi},'same'); % 1xtime complex
    ppcMat(fi,fi,:) = tmp(dvec); %subsample every 10th value from tmp
    tmp = conv(pulLFP,wavs{fi},'same');
    pulMat(fi,fi,:) = tmp(dvec);
    tmp = conv(vcLFP,wavs{fi},'same');
    vcMat(fi,fi,:) = tmp(dvec);
    for fj = fi+1:numFreqs
        tmp  = conv(ppcLFP,wavs{fj},'same');
        tmp2 = conv(abs(tmp),wavs{fi},'same');
        ppcMat(fi,fj,:) = tmp2(dvec);
        tmp  = conv(pulLFP,wavs{fj},'same');
        tmp2 = conv(abs(tmp),wavs{fi},'same');
        pulMat(fi,fj,:) = tmp2(dvec);
        tmp  = conv(vcLFP,wavs{fj},'same');
        tmp2 = conv(abs(tmp),wavs{fi},'same');
        vcMat(fi,fj,:) = tmp2(dvec);
    end
end


% within region
% compute coherence
for fi = 1:numFreqs
     display(num2str(fi))
    for fj = fi+1:numFreqs
        % PPC
        Slh = mean(squeeze(ppcMat(fi,fi,:)).*conj(squeeze(ppcMat(fi,fj,:))));
        Sll = mean(squeeze(ppcMat(fi,fi,:)).*conj(squeeze(ppcMat(fi,fi,:))));
        Shh = mean(squeeze(ppcMat(fi,fj,:)).*conj(squeeze(ppcMat(fi,fj,:))));
        ppc_Cy(fi,fj)  = Slh/sqrt(Sll*Shh);
       
        r = corrcoef(squeeze(abs(ppcMat(fi,fi,:))),squeeze(abs(ppcMat(fj,fj,:))));
        ppc_r(fi,fj) = r(1,2);
       
        % LP/Pulvinar
        Slh = mean(squeeze(pulMat(fi,fi,:)).*conj(squeeze(pulMat(fi,fj,:))));
        Sll = mean(squeeze(pulMat(fi,fi,:)).*conj(squeeze(pulMat(fi,fi,:))));
        Shh = mean(squeeze(pulMat(fi,fj,:)).*conj(squeeze(pulMat(fi,fj,:))));
        pul_Cy(fi,fj)  = Slh/sqrt(Sll*Shh);
       
        r = corrcoef(squeeze(abs(pulMat(fi,fi,:))),squeeze(abs(pulMat(fj,fj,:))));
        pul_r(fi,fj) = r(1,2);
       
        % visual cortex
        Slh = mean(squeeze(vcMat(fi,fi,:)).*conj(squeeze(vcMat(fi,fj,:))));
        Sll = mean(squeeze(vcMat(fi,fi,:)).*conj(squeeze(vcMat(fi,fi,:))));
        Shh = mean(squeeze(vcMat(fi,fj,:)).*conj(squeeze(vcMat(fi,fj,:))));
        vc_Cy(fi,fj)  = Slh/sqrt(Sll*Shh);
       
        r = corrcoef(squeeze(abs(vcMat(fi,fi,:))),squeeze(abs(vcMat(fj,fj,:))));
        vc_r(fi,fj) = r(1,2);
    end
end

save([GroupAnalysisDir 'PAC_local.mat'],'pul_Cy','pul_r','ppc_Cy','ppc_r','vc_Cy','vc_r', '-v7.3');

%% plot PAC and significance
fois = [2 4 8 16 32 64 128]; % find the index of each foi that is closest to each fois, i.e. abs(foi-fois)~0
for fi = 1:numel(fois)
    [bi,bb] = sort(abs(foi-fois(fi)));
    tickLoc(fi) = bb(1);
end


% Cross frequency coupling
 
% PPC
tickLabel = {'2','4','8','16','32','64','128'};

fig = figure();
% LP/Pulvinar
subplot(3,3,1)
imagesc(1:numel(foi),1:numel(foi),flipud(rot90((abs(pul_Cy)))));
set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',tickLabel)
xlabel('LP/Pulvinar phase frequency')
ylabel('LP/Pulvinar amplitude frequency')
colormap(awesomeMap)
colorbar;
%caxis([0 0.5])
freezeColors
title('Phase amplitude coupling')
 
subplot(3,3,2)
imagesc(1:numel(foi),1:numel(foi),flipud(rot90(pul_r)));
set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',tickLabel)
xlabel('LP/Pulvinar phase frequency');
ylabel('LP/Pulvinar amplitude frequency');
colormap(bry_map)
colorbar;
%caxis([-0.2 0.2])
title('Amplitude correlation')


%PPC
subplot(3,2,3)
imagesc(1:numel(foi),1:numel(foi),flipud(rot90((abs(ppc_Cy)))));
set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',tickLabel)
xlabel('PPC phase frequency')
ylabel('PPC amplitude frequency')
colormap(awesomeMap)
colorbar;
%caxis([0 0.7])
freezeColors
title('Phase amplitude coupling')
 
subplot(3,2,4)
imagesc(1:numel(foi),1:numel(foi),flipud(rot90(ppc_r)));
set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',tickLabel)
xlabel('PPC phase frequency')
ylabel('PPC amplitude frequency')
colormap(bry_map)
colorbar;
%caxis([-0.4 0.4])
title('Amplitude correlation')
 
% Visual cortex
subplot(3,2,5)
imagesc(1:numel(foi),1:numel(foi),flipud(rot90((abs(vc_Cy)))));
set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',tickLabel)
xlabel('Visual cortex phase frequency');
ylabel('Visual cortex amplitude frequency');
colormap(awesomeMap);
colorbar;
%caxis([0 0.7])
freezeColors
title('Phase amplitude coupling')
 
subplot(3,2,6)
imagesc(1:numel(foi),1:numel(foi),flipud(rot90(vc_r)));
set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',tickLabel)
xlabel('Visual cortex phase frequency');
ylabel('Visual cortex amplitude frequency');
colormap(bry_map);
colorbar;
%caxis([-0.4 0.4])
title('Amplitude correlation')

save([GroupAnalysisDir, 'PAC_2-7min_a4b14.mat'],'foi','ppc_Cy','ppc_r','vc_Cy','vc_r','pul_Cy','pul_r','-v7.3');

savefig(fig, [GroupAnalysisDir, 'PAC_2-7min_a4b14.fig']);
saveas(fig,[GroupAnalysisDir, 'PAC_2-7min_a4b14.png']);

%% plot PAC only for frequency range of interest
% rescale
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-150)/2 (screensize(4)-150)/4]); %(x,y,width,height) screensize(3)-100

%xLim = [1,75];yLim = [75,150]; %1-16Hz numFreq = 150, half is 16Hz, 75th
xFOI = [2,8]; yFOI = [8,128];
xLim = [find(foi>=xFOI(1),1), find(foi>=xFOI(2),1)];
yLim = [find(foi>=yFOI(1),1), length(foi)];

subplot(1,3,1)
imagesc(1:numel(foi),1:numel(foi),flipud(rot90((abs(pul_Cy)))));
set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',tickLabel)
xlabel('Phase frequency [Hz]');
ylabel('Amplitude frequency [Hz]');
xlim(xLim);ylim(yLim);
%colormap(awesomeMap)
colorbar;
caxis([0 0.8])
%freezeColors
title('LPl')

%PPC
subplot(1,3,2)
imagesc(1:numel(foi),1:numel(foi),flipud(rot90((abs(ppc_Cy)))));
set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',tickLabel)
xlabel('Phase frequency [Hz]');
ylabel('Amplitude frequency [Hz]');
xlim(xLim);ylim(yLim);
% colormap(awesomeMap)
colorbar;
caxis([0 0.4])
%freezeColors
title('PPC')

% Visual cortex
subplot(1,3,3)
imagesc(1:numel(foi),1:numel(foi),flipud(rot90((abs(vc_Cy)))));
set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',tickLabel)
xlabel('Phase frequency [Hz]');
ylabel('Amplitude frequency [Hz]');
xlim(xLim);ylim(yLim);
% colormap(awesomeMap);
colorbar;
caxis([0 0.4])
% freezeColors
title('V1')
colormap(awesomeMap)

% subplot(1,4,4)
% caxis([0 0.6])
% colorbar;

savefig(fig, [GroupAnalysisDir, 'PAC_2-7min_a4b14_' num2str(xFOI(1)) '-' num2str(xFOI(2)) 'Hz.fig']);
saveas(fig,[GroupAnalysisDir, 'PAC_2-7min_a4b14_' num2str(xFOI(1)) '-' num2str(xFOI(2)) 'Hz.png']);
    end
end
x = toc;
fprintf('time required =%f sec\n', x);