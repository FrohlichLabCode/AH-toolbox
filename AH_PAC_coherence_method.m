% This code is adopted from Iain's working code, the result doesn't match
% what we see in the PAC_PLV_method, which is the conventional way of
% calculating PAC and also the way we used in our paper.
% So be careful when using this code. 


clear

cluster = 0;
tic

%% load data and functions

if cluster == 0
    addpath('E:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys'); % all the functions needed
    %addpath('J:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\mvgc_v1.0') % for Granger Causality
    ephysDir     = ['E:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\Analyzed\0153_AttentionTask3_11\'];     
    saveRootPath = ['E:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\Analyzed\0153_AttentionTask3_11\PAC\'];
    mkdir(saveRootPath);  
    cd(ephysDir); % change directory to save data

    % load lfpData
    load([ephysDir 'lfpPCA.mat']);
        
elseif cluster == 1
    % % Required for use of KilLDevil and parfor
    % ClusterInfo.setQueueName('week')
    % ClusterInfo.setUserDefinedOptions('-M72')
    % ClusterInfo.setUserDefinedOptions('-R mem128')
    % % add paths for cluster computation
    addpath(genpath('/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
    pathDir      = [ '/pine/scr/a/n/angelvv/5CSRTT/Analyzed/0153_AttentionTask3_14']; % CHANGE FOR KILLDEVIL VS LONGLEAF
    saveRootPath = [ '/pine/scr/a/n/angelvv/5CSRTT/Analyzed/0153_AttentionTask3_14/PAC'];
    mkdir(saveRootPath)
    %code for initialising parallel computing
    
    numCore = 24; % USR DEFINE
    parpool('local',numCore);
    fileInfo     = dir([pathDir '0*']);

    ephysDir = '/pine/scr/h/w/angelvv/5CSRTT/Analyzed/0153_AttentionTask3_11/';
    cd(ephysDir); % change directory to save data
    addpath('/nas/longleaf/home/angelvv/Code/'); % all the functions needed
    addpath('/nas/longleaf/home/angelvv/Code/mvgc_v1.0/')
    
    % load lfpData
    load([ephysDir '/lfpPCA.mat']);
    
end

%% 
% analyze using lfpPCA, 1 row for each region

lowFreq      = 2;  % wavelet is significantly slow when freq is <1
highFreq     = 80;
numFreqs     = 50;
lfpFs        = 1000;
foi          = logspace(log10(lowFreq),log10(highFreq),numFreqs); %generates a row vector logarithmically spaced points between decades 10^a and 10^b

lfpMat = lfpPCA(:,2*60*lfpFs:7*60*lfpFs); %example

% downsample to 200Hz
newFs   = 200;
[nr,dr] = rat(lfpFs/newFs); %returns two arrays, N and D, such that N./D ~ X
downMat = resample(lfpMat',dr,nr)'; %resamples matrix at p/q times the original sample rate. If x is a matrix, then resample treats each column of x as an independent channe

ppcChans = 1;
vcChans  = 3;
pulChans = 2;

ppcLFP = nanmedian(downMat(ppcChans,:),1); % here only 1 PCA channel but can also use median of many channels
vcLFP  = nanmedian(downMat(vcChans,:),1);
pulLFP = nanmedian(downMat(pulChans,:),1);

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

fois = [2 4 8 16 32 64]; % find the index of each foi that is closest to each fois, i.e. abs(foi-fois)~0
for fi = 1:numel(fois)
    [bi,bb] = sort(abs(foi-fois(fi)));
    tickLoc(fi) = bb(1);
end
 
% Cross frequency coupling
 
% PPC
tickLabel = {'2','4','8','16','32','64'};

fig = figure();
subplot(3,2,1)
imagesc(1:numel(foi),1:numel(foi),flipud(rot90((abs(ppc_Cy)))));
set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',tickLabel)
xlabel('PPC phase frequency')
ylabel('PPC amplitude frequency')
colormap(awesomeMap)
colorbar;
%caxis([0 0.7])
freezeColors
title('Phase amplitude coupling')
 
subplot(3,2,2)
imagesc(1:numel(foi),1:numel(foi),flipud(rot90(ppc_r)));
set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',tickLabel)
xlabel('PPC phase frequency')
ylabel('PPC amplitude frequency')
colormap(bry_map)
colorbar;
%caxis([-0.4 0.4])
title('Amplitude correlation')
 
% LP/Pulvinar
subplot(3,2,3)
imagesc(1:numel(foi),1:numel(foi),flipud(rot90((abs(pul_Cy)))));
set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',tickLabel)
xlabel('LP/Pulvinar phase frequency')
ylabel('LP/Pulvinar amplitude frequency')
colormap(awesomeMap)
colorbar;
%caxis([0 0.5])
freezeColors
title('Phase amplitude coupling')
 
subplot(3,2,4)
imagesc(1:numel(foi),1:numel(foi),flipud(rot90(pul_r)));
set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',tickLabel)
xlabel('LP/Pulvinar phase frequency');
ylabel('LP/Pulvinar amplitude frequency');
colormap(bry_map)
colorbar;
%caxis([-0.2 0.2])
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

cd(saveRootPath);
savefig(fig, ['PAC_2-7min.fig']);
save('PAC data_2-7min.mat','ppc_Cy','ppc_r','vc_Cy','vc_r','pul_Cy','pul_r','-v7.3');

x = toc;
fprintf('time required =%f sec\n', x);