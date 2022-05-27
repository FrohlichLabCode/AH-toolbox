
%% 1.1 Calculate phase-amplitude-coupling based PLV method for entire recording

globalPath = 'J:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\Analyzed\';
recs       = dir([globalPath '0153_AttentionTask3*\1-80Hzlog\funcConPCA_region1-2.mat']);

for irec = 1:numel(recs)
    display(['rec ' num2str(irec) '\' num2str(numel(recs))])
%     recName = recs(irec).folder(end-1:end);
%     
%     if exist([globalPath '\PAC\' recName '_PAC_spectrum.mat'],'file'); continue; end
%     % load in meta data and digital triggers
%     load([globalPath recName '_MetaData'])
%     load([globalPath recName '_Digital'])

    load([recs(irec).folder '\' recs(irec).name]);%(output from AH_PCA_functional_connectivity.m)
    vars = {'funcCon', 'GC_XtoY_All', 'GC_YtoX_All', 'iCoherenceAll','imagZAll','plvAll','psiAll','psiNormAll','xSpecAll','ySpecAll'};
    clear(vars{:}); % free up memory 
    
%     %%% find trials and event onsets
%     StimOnset            = Digital.dig3On/Digital.Fs;
%     SelectedTrialCorrect = MetaData.SelectedTrialCorrect;
%     StimOnsetCorrect     = StimOnset( SelectedTrialCorrect);
%     chans2analyse        = MetaData.SelectedChn;
    
%     % load in LFP
%     load([globalPath recName '_LFP'])
%     lfpMat = LFP.LFPdata(chans2analyse,:);
%     lfpFs  = LFP.Fs/LFP.downsample;
    
    % downsample to 200Hz
    lfpMat = lfpPCA;
    newFs   = 200;
    [nr,dr] =  rat(lfpFs/newFs);
    for ichan = 1:size(lfpMat,1) 
        downMat(ichan,:) = resample(lfpMat(ichan,:)',dr,nr)';  % double-checked
    end
    
    lowFreqs  = 1:0.25:10;
    highFreqs = 15:80;
    lowWav    = is_makeWavelet(lowFreqs,newFs);
    highWav   = is_makeWavelet(highFreqs,newFs);
    
    for ichan = 1:size(downMat,1)
        testData = downMat(ichan,:); %1 x Timepoints
        lowMat   = nan(numel(lowFreqs),numel(testData));
        highMat  = nan(numel(highFreqs),numel(lowFreqs),numel(testData));
        
        for lf = 1:numel(lowFreqs)
            lowMat(lf,:) = conv(testData,lowWav{lf},'same');
        end
        
        for hf = 1:numel(highFreqs)
            highTmp = conv(testData,highWav{hf},'same');
            highAmp = abs(highTmp);
            for lf = 1:numel(lowFreqs)
                highMat(hf,lf,:) = conv(highAmp,lowWav{lf},'same');
            end
        end
        
        lowAng  = angle(lowMat);
        highAng = angle(highMat);
        
        for lf = 1:numel(lowFreqs)
            for hf = 1:numel(highFreqs)
                ang = lowAng(lf,:)' - squeeze(highAng(hf,lf,:));
                plv(ichan,hf,lf) = abs(nanmean(exp(1i*ang)));
            end
        end
    end
    
    save([recs(irec).folder '\PAC\PAC_spectrum.mat'],'plv','lowFreqs','highFreqs');
    vars = {'downMat', 'lfpMat'};
    clear(vars{:});
end

%% 1.2 plot each session's PAC
globalPath = 'J:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\Analyzed\';
recs       = dir([globalPath '0153_AttentionTask3*\PAC\PAC_spectrum.mat']);
PAC_spectrum_all    = []; % size(PAC_spectrum_all) = [3 4 66 37]

for irec = 1:numel(recs)
    display(['rec ' num2str(irec) '\' num2str(numel(recs))])
    
    load([recs(irec).folder '\' recs(irec).name]);
    PAC_spectrum_all(:,irec,:,:) = plv;
    
    doPlot = 0;
    if doPlot == 1
    screensize = get( groot, 'Screensize');
    fig = figure('Position',[50 50 screensize(3)-100 screensize(4)/2.5]);
    subplot(1,3,1)
    imagesc(lowFreqs,highFreqs,squeeze(plv(1,:,:))) % 1st channel is PPC
    set(gca,'Ydir','normal','TickDir','out')
    colormap(awesomeMap)
    c = colorbar; ylabel(c,'PAC (PLV)');
    caxis([0 0.7])
    xlabel('Phase frequency')
    ylabel('Amplitude frequency')
    title('PPC phase amplitude coupling')

    subplot(1,3,2)
    imagesc(lowFreqs,highFreqs,squeeze(plv(2,:,:)))
    set(gca,'Ydir','normal','TickDir','out')
    colormap(awesomeMap)
    c = colorbar; ylabel(c,'PAC (PLV)');
    caxis([0 0.7])
    xlabel('Phase frequency')
    ylabel('Amplitude frequency')
    title('LPl phase amplitude coupling')

    subplot(1,3,3)
    imagesc(lowFreqs,highFreqs,squeeze(plv(3,:,:)))
    set(gca,'Ydir','normal','TickDir','out')
    colormap(awesomeMap)
    c = colorbar; ylabel(c,'PAC (PLV)');
    caxis([0 0.7])
    xlabel('Phase frequency')
    ylabel('Amplitude frequency')
    title('VC phase amplitude coupling')

    saveDir = [recs(irec).folder '\'];
    savefig(fig,[saveDir 'low_high_PAC'])
    saveas(fig, [saveDir 'low_high_PAC.png']);
    end
end

save([globalPath 'PAC_spectrum_all.mat'],'PAC_spectrum_all','lowFreqs','highFreqs');

% plot the average PAC
% normalize PAC across sessions
for iRegion = 1:size(PAC_spectrum_all,1)
    temp = squeeze(PAC_spectrum_all(iRegion,:,:,:));
    session_mean = nanmean(nanmean(temp,2),3);
    all_mean = nanmean(session_mean);
    for iSession = 1:size(PAC_spectrum_all,2)
        PAC_normed(iRegion,iSession,:,:) = squeeze(temp(iSession,:,:))./session_mean(iSession).*all_mean;
    end
end

PAC_spec_median = squeeze(nanmedian(PAC_normed,2));

doPlot = 1;
if doPlot == 1
    screensize = get( groot, 'Screensize');
    fig = figure('Position',[50 50 screensize(3)-100 screensize(4)/2.5]);
    subplot(1,3,1)
    imagesc(lowFreqs,highFreqs,squeeze(PAC_spec_median(1,:,:))) % 1st channel is PPC
    set(gca,'Ydir','normal','TickDir','out')
    colormap(awesomeMap)
    c = colorbar; ylabel(c,'PAC (PLV)');
    caxis([0 0.7])
    xlabel('Phase frequency')
    ylabel('Amplitude frequency')
    title('PPC phase amplitude coupling')

    subplot(1,3,2)
    imagesc(lowFreqs,highFreqs,squeeze(PAC_spec_median(2,:,:)))
    set(gca,'Ydir','normal','TickDir','out')
    colormap(awesomeMap)
    c = colorbar; ylabel(c,'PAC (PLV)');
    caxis([0 0.7])
    xlabel('Phase frequency')
    ylabel('Amplitude frequency')
    title('LPl phase amplitude coupling')

    subplot(1,3,3)
    imagesc(lowFreqs,highFreqs,squeeze(PAC_spec_median(3,:,:)))
    set(gca,'Ydir','normal','TickDir','out')
    colormap(awesomeMap)
    c = colorbar; ylabel(c,'PAC (PLV)');
    caxis([0 0.7])
    xlabel('Phase frequency')
    ylabel('Amplitude frequency')
    title('VC phase amplitude coupling')

    savefig(fig,[globalPath 'PAC_spectrum_median'])
    saveas(fig, [globalPath 'PAC_spectrum_median.png']);
end

%% 2.1 based on PAC plot, sum over frequences of interest and put all sessions together 
% get sum within the range for power and median for plv and coherence
lf_range = [4,8];  % theta
hf_range = [40,80];% gamma
timeWin  = [-4,2];

PPC_pow_lf = []; PPC_pow_hf = [];
LPl_pow_lf = []; LPl_pow_hf = [];
VC_pow_lf  = []; VC_pow_hf  = [];
PPC_LPl_plv_lf = []; PPC_LPl_plv_hf = [];
PPC_LPl_coh_lf = []; PPC_LPl_coh_hf = [];
PPC_VC_plv_lf  = []; PPC_VC_plv_hf  = [];
PPC_VC_coh_lf  = []; PPC_VC_coh_hf  = [];
LPl_VC_plv_lf  = []; LPl_VC_plv_hf  = [];
LPl_VC_coh_lf  = []; LPl_VC_coh_hf  = [];


lfMask = foi>lf_range(1) & foi<=lf_range(2);
hfMask = foi>hf_range(1) & foi<=hf_range(2);

%1-2
globalPath = 'J:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\Analyzed\';
recs       = dir([globalPath '0153_AttentionTask3*\1-80Hzlog\funcConPCA_region1-2.mat']);
for irec = 1:numel(recs)
    display(['region1-2 rec ' num2str(irec) '\' num2str(numel(recs))])
    
    load([recs(irec).folder '\' recs(irec).name]);
    PPC_pow_lf = [PPC_pow_lf; nansum(squeeze(funcCon.xspec(lfMask,:)),1)];
    PPC_pow_hf = [PPC_pow_hf; nansum(squeeze(funcCon.xspec(hfMask,:)),1)];
    LPl_pow_lf = [LPl_pow_lf; nansum(squeeze(funcCon.yspec(lfMask,:)),1)];
    LPl_pow_hf = [LPl_pow_hf; nansum(squeeze(funcCon.yspec(hfMask,:)),1)];
    PPC_LPl_plv_lf = [PPC_LPl_plv_lf; nanmean(squeeze(funcCon.plv(lfMask,:)),1)];
    PPC_LPl_plv_hf = [PPC_LPl_plv_hf; nanmean(squeeze(funcCon.plv(hfMask,:)),1)];
    PPC_LPl_coh_lf = [PPC_LPl_coh_lf; nanmean(squeeze(funcCon.coherence(lfMask,:)),1)];
    PPC_LPl_coh_hf = [PPC_LPl_coh_hf; nanmean(squeeze(funcCon.coherence(hfMask,:)),1)];
end

%1-3
recs       = dir([globalPath '0153_AttentionTask3*\1-80Hzlog\funcConPCA_region1-3.mat']);
for irec = 1:numel(recs)
    display(['region1-3 rec ' num2str(irec) '\' num2str(numel(recs))])
    
    load([recs(irec).folder '\' recs(irec).name]);
    VC_pow_lf = [VC_pow_lf; nansum(squeeze(funcCon.yspec(lfMask,:)),1)];
    VC_pow_hf = [VC_pow_hf; nansum(squeeze(funcCon.yspec(hfMask,:)),1)];
    PPC_VC_plv_lf = [PPC_VC_plv_lf; nanmean(squeeze(funcCon.plv(lfMask,:)),1)];
    PPC_VC_plv_hf = [PPC_VC_plv_hf; nanmean(squeeze(funcCon.plv(hfMask,:)),1)];
    PPC_VC_coh_lf = [PPC_VC_coh_lf; nanmean(squeeze(funcCon.coherence(lfMask,:)),1)];
    PPC_VC_coh_hf = [PPC_VC_coh_hf; nanmean(squeeze(funcCon.coherence(hfMask,:)),1)];
end


%2-3
recs       = dir([globalPath '0153_AttentionTask3*\1-80Hzlog\funcConPCA_region2-3.mat']);
for irec = 1:numel(recs)
    display(['region2-3 rec ' num2str(irec) '\' num2str(numel(recs))])
    
    load([recs(irec).folder '\' recs(irec).name]);
    LPl_VC_plv_lf = [LPl_VC_plv_lf; nanmean(squeeze(funcCon.plv(lfMask,:)),1)];
    LPl_VC_plv_hf = [LPl_VC_plv_hf; nanmean(squeeze(funcCon.plv(hfMask,:)),1)];
    LPl_VC_coh_lf = [LPl_VC_coh_lf; nanmean(squeeze(funcCon.coherence(lfMask,:)),1)];
    LPl_VC_coh_hf = [LPl_VC_coh_hf; nanmean(squeeze(funcCon.coherence(hfMask,:)),1)];
end

%% 3.1 stats, normalize by session total power and get median and std
PPC_lf_stats = medianStd(PPC_pow_lf);
PPC_hf_stats = medianStd(PPC_pow_hf);
LPl_lf_stats = medianStd(LPl_pow_lf);
LPl_hf_stats = medianStd(LPl_pow_hf);
VC_lf_stats  = medianStd(VC_pow_lf);
VC_hf_stats  = medianStd(VC_pow_hf);

PPC_LPl_plv_lf_stats = medianStd(PPC_LPl_plv_lf);
PPC_LPl_plv_hf_stats = medianStd(PPC_LPl_plv_hf);
PPC_VC_plv_lf_stats  = medianStd(PPC_VC_plv_lf);
PPC_VC_plv_hf_stats  = medianStd(PPC_VC_plv_hf);
LPl_VC_plv_lf_stats  = medianStd(LPl_VC_plv_lf);
LPl_VC_plv_hf_stats  = medianStd(LPl_VC_plv_hf);

PPC_LPl_coh_lf_stats = medianStd(PPC_LPl_coh_lf);
PPC_LPl_coh_hf_stats = medianStd(PPC_LPl_coh_hf);
PPC_VC_coh_lf_stats  = medianStd(PPC_VC_coh_lf);
PPC_VC_coh_hf_stats  = medianStd(PPC_VC_coh_hf);
LPl_VC_coh_lf_stats  = medianStd(LPl_VC_coh_lf);
LPl_VC_coh_hf_stats  = medianStd(LPl_VC_coh_hf);
save([globalPath '\FreqBandMean\FreqBandMean_4-8Hz_40-80Hz.mat'],'LPl*','PPC*','VC*');
   
%% 4 plot spectrum for lf band and hf band
% 4.1 power spectrum for each region
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 screensize(3)-100 (screensize(4)-150)/2]);
subplot(1,3,3)
x  = funcCon.tvec';
y  = LPl_lf_stats(1,:)';
dy = LPl_lf_stats(2,:)';  % made-up error values
low = y-dy; low(low<0) = 0.1; %Power has to >0 to take log (convert to db)
h  = fill([x;flipud(x)],[pow2db(low);pow2db(flipud(y+dy))],'b','linestyle','none');
set(h,'facealpha',.5); % make translucent
l1 = line(x,pow2db(y),'Color','blue');

hold on
y  = LPl_hf_stats(1,:)';
dy = LPl_hf_stats(2,:)';  % made-up error values
low = y-dy; low(low<0) = 0.1;
h  = fill([x;flipud(x)],[pow2db(low);pow2db(flipud(y+dy))],'r','linestyle','none');
set(h,'facealpha',.5); % make translucent
l2 = line(x,pow2db(y),'Color','red');
legend([l1 l2],'theta','gamma','Location','east');
%ylim([65,75])
xlabel('Time to touch [sec]')
ylabel('LPl power [dB(uV^2)]')

subplot(1,3,1)
x  = funcCon.tvec';
y  = PPC_lf_stats(1,:)';
dy = PPC_lf_stats(2,:)';  % made-up error values
low = y-dy; low(low<0) = 0.1; %Power has to >0 to take log (convert to db)
h  = fill([x;flipud(x)],[pow2db(low);pow2db(flipud(y+dy))],'b','linestyle','none');
set(h,'facealpha',.5); % make translucent
l1 = line(x,pow2db(y),'Color','blue');

hold on
y  = PPC_hf_stats(1,:)';
dy = PPC_hf_stats(2,:)';  % made-up error values
low = y-dy; low(low<0) = 0.1;
h  = fill([x;flipud(x)],[pow2db(low);pow2db(flipud(y+dy))],'r','linestyle','none');
set(h,'facealpha',.5); % make translucent
l2 = line(x,pow2db(y),'Color','red');
%legend([l1 l2],'theta','gamma')
ylim([20,70])
xlabel('Time to touch [sec]')
ylabel('PPC power [dB(uV^2)]')


subplot(1,3,2)
x  = funcCon.tvec';
y  = VC_lf_stats(1,:)';
dy = VC_lf_stats(2,:)';  % made-up error values
low = y-dy; low(low<0) = 0.1;
h  = fill([x;flipud(x)],[pow2db(low);pow2db(flipud(y+dy))],'b','linestyle','none');
set(h,'facealpha',.5); % make translucent
l1 = line(x,pow2db(y),'Color','blue');

hold on
y  = VC_hf_stats(1,:)';
dy = VC_hf_stats(2,:)';  % made-up error values
low = y-dy; low(low<0) = 0.1; 
h  = fill([x;flipud(x)],[pow2db(low);pow2db(flipud(y+dy))],'r','linestyle','none');
set(h,'facealpha',.5); % make translucent
l2 = line(x,pow2db(y),'Color','red');
%legend([l1 l2],'theta','gamma')
%ylim([65,75])
xlabel('Time to touch [sec]')
ylabel('VC power [dB(uV^2)]')

savefig(fig,[globalPath '\FreqBandMean\FreqBandMean_4-8Hz_40-80Hz_pow_auto.fig'],'compact');
saveas(fig,[globalPath '\FreqBandMean\FreqBandMean_4-8Hz_40-80Hz_pow_auto.png']);

%% 4.2 plv over time for each region pair
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 screensize(3)-100 (screensize(4)-150)/2.5]);
subplot(1,3,1)
x  = funcCon.tvec;
y  = PPC_LPl_plv_lf_stats;
H1 = shadedErrorBar(x,y(1,:),y(2,:),'-b',1); 
hold on

y  = PPC_LPl_plv_hf_stats;
H2 = shadedErrorBar(x,y(1,:),y(2,:),'-r',1); 
xlabel('Time to touch [sec]')
ylabel('PPC/LPl PLV')
legend([H1.mainLine H2.mainLine],'theta','gamma')
ylim([0,1])
hold off

subplot(1,3,2)
y  = PPC_VC_plv_lf_stats;
H1 = shadedErrorBar(x,y(1,:),y(2,:),'-b',1); 
hold on

y  = PPC_VC_plv_hf_stats;
H2 = shadedErrorBar(x,y(1,:),y(2,:),'-r',1); 
xlabel('Time to touch [sec]')
ylabel('PPC/VC PLV')
%legend([H1.mainLine H2.mainLine],'theta','gamma')
ylim([0,1])
hold off

subplot(1,3,3)
y  = LPl_VC_plv_lf_stats;
H1 = shadedErrorBar(x,y(1,:),y(2,:),'-b',1); 
hold on

y  = LPl_VC_plv_hf_stats;
H2 = shadedErrorBar(x,y(1,:),y(2,:),'-r',1); 
xlabel('Time to touch [sec]')
ylabel('LPl/VC PLV')
%legend([H1.mainLine H2.mainLine],'theta','gamma')
ylim([0,1])
hold off

savefig(fig,[globalPath '\FreqBandMean\FreqBandMean_4-8Hz_40-80Hz_plv_same.fig'],'compact');
saveas(fig,[globalPath '\FreqBandMean\FreqBandMean_4-8Hz_40-80Hz_plv_same.png']);


%% 4.3 Coherence over time for each region pair
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 screensize(3)-100 (screensize(4)-150)/2.5]);
subplot(1,3,1)
x  = funcCon.tvec;
y  = PPC_LPl_coh_lf_stats;
H1 = shadedErrorBar(x,y(1,:),y(2,:),'-b',1); 
hold on

y  = PPC_LPl_coh_hf_stats;
H2 = shadedErrorBar(x,y(1,:),y(2,:),'-r',1); 
xlabel('Time to touch [sec]')
ylabel('PPC/LPl coherence')
legend([H1.mainLine H2.mainLine],'theta','gamma')
ylim([0,1])
hold off

subplot(1,3,2)
y  = PPC_VC_coh_lf_stats;
H1 = shadedErrorBar(x,y(1,:),y(2,:),'-b',1); 
hold on

y  = PPC_VC_coh_hf_stats;
H2 = shadedErrorBar(x,y(1,:),y(2,:),'-r',1); 
xlabel('Time to touch [sec]')
ylabel('PPC/VC coherence')
%legend([H1.mainLine H2.mainLine],'theta','gamma')
ylim([0,1])
hold off

subplot(1,3,3)
y  = LPl_VC_coh_lf_stats;
H1 = shadedErrorBar(x,y(1,:),y(2,:),'-b',1); 
hold on

y  = LPl_VC_coh_hf_stats;
H2 = shadedErrorBar(x,y(1,:),y(2,:),'-r',1); 
xlabel('Time to touch [sec]')
ylabel('LPl/VC coherence')
%legend([H1.mainLine H2.mainLine],'theta','gamma')
ylim([0,1])
hold off

savefig(fig,[globalPath '\FreqBandMean\FreqBandMean_4-8Hz_40-80Hz_coh_same.fig'],'compact');
saveas(fig,[globalPath '\FreqBandMean\FreqBandMean_4-8Hz_40-80Hz_coh_same.png']);

%% 4.4 stats on theta, gamma
% PPC_LPl_plv = [];
% PPC_LPl_coh = [];
% PPC_VC_plv = [];
% PPC_VC_coh = [];
% LPl_VC_plv = [];
% LPl_VC_coh = [];



%% 6 Time resolved PAC % this takes ~20min on GPU computer
globalPath = 'J:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\Analyzed\';
recs       = dir([globalPath '0153_AttentionTask3*\1-80Hzlog\lfpPCA.mat']);
lf_range = [4,7];
hf_range = [45,80];
tPAC_all_session = []; % store PAC mean (PPC1,LPl1,VC1,PPC2,LPl2...)
for irec = 1:numel(recs)
    display(['rec ' num2str(irec) '\' num2str(numel(recs))])
    load([recs(irec).folder '\' recs(irec).name]);
    saveRootPath = recs(irec).folder;
    tPAC = PAC_time_resolved(lfpPCA,lfpFs,event,twin,saveRootPath);
    lf_mask = tPAC.lowFreqs>lf_range(1) & tPAC.lowFreqs<lf_range(2);
    hf_mask = tPAC.highFreqs>hf_range(1) & tPAC.highFreqs<hf_range(2);
    for iregion = 1:size(lfpPCA,1)
        tPAC_all_session = [tPAC_all_session; squeeze(nanmean(nanmean(tPAC.plv{iregion}(hf_mask, lf_mask,:),1),2))'];
    end
end

tPAC_stats = struct;
nregion = size(lfpPCA,1);
nrec = size(recs,1);
tPAC_stats.PPC = medianStd(tPAC_all_session(1:nregion:nregion*nrec,:));
tPAC_stats.LPl = medianStd(tPAC_all_session(2:nregion:nregion*nrec,:));
tPAC_stats.VC  = medianStd(tPAC_all_session(3:nregion:nregion*nrec,:));
save([globalPath 'tPAC_4-7Hz_45-80Hz.mat'],'tPAC_all_session','tPAC_stats');

% plot tPAC for each region
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 screensize(3)-100 (screensize(4)-150)/2.5]);
subplot(1,3,1)
x  = linspace(twin(1),twin(2),size(tPAC_all_session,2));
y  = tPAC_stats.PPC;
shadedErrorBar(x,y(1,:),y(2,:),'-k',1); 
ylim([0.1,0.7])
xlabel('Time to touch [sec]')
ylabel('PPC theta/gamma PAC')

subplot(1,3,2)
x  = linspace(twin(1),twin(2),size(tPAC_all_session,2));
y  = tPAC_stats.VC;
shadedErrorBar(x,y(1,:),y(2,:),'-k',1); 
ylim([0.1,0.7])
xlabel('Time to touch [sec]')
ylabel('VC theta/gamma PAC')

subplot(1,3,3)
x  = linspace(twin(1),twin(2),size(tPAC_all_session,2));
y  = tPAC_stats.LPl;
shadedErrorBar(x,y(1,:),y(2,:),'-k',1); 
ylim([0.1,0.7])
xlabel('Time to touch [sec]')
ylabel('LPl theta/gamma PAC')

savefig(fig,[globalPath '\tPAC_4-7Hz_45-80Hz.fig'],'compact');
saveas(fig,[globalPath '\tPAC_4-7Hz_45-80Hz.png']);

% plot tPAC for 3 regions together
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-100)/5 (screensize(4)-150)/2.5]);
x  = linspace(twin(1),twin(2),size(tPAC_all_session,2));
y  = tPAC_stats.PPC;
H1 = shadedErrorBar(x,y(1,:),y(2,:),'-r',1); 
hold on
y  = tPAC_stats.VC;
H2 = shadedErrorBar(x,y(1,:),y(2,:),'-b',1); 
hold on
y  = tPAC_stats.LPl;
H3 = shadedErrorBar(x,y(1,:),y(2,:),'-m',1);
vline(0,'--k');
legend([H1.mainLine H2.mainLine H3.mainLine],'PPC','V1','LPl')%,'Location','east');
ylim([0.1,0.75])
xlabel('Time to touch [sec]')
ylabel('Theta/gamma PAC')
hold off
savefig(fig,[globalPath '\tPAC_4-7Hz_45-80Hz_3in1.fig'],'compact');
saveas(fig,[globalPath '\tPAC_4-7Hz_45-80Hz_3in1.png']);