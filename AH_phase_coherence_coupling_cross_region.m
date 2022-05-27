% recPath = 'D:\ferret_working_dir\0153\rawData\toAnalyze\0153_AttentionTask2Ephys_01_20171020_171020_141018\';
% load([recPath 'lfp\lfpMat'])
load('J:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\Analyzed\0153_AttentionTask3_13\1-80Hzlog\lfpPCA.mat');
lfpMat = lfpPCA;
addpath('J:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys'); % all the functions needed
addpath(genpath('J:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\mvgc_v1.0')); % use ptic, ptoc

saveRootDir = 'J:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\Analyzed\0153_AttentionTask3_13\1-80Hzlog\';
%% Plot LFP
% figure()
% for ichan = 1:16
%     load([recPath 'spikes\spk_' num2str(ichan+64)]);
%  
%     subplot(3,6,ichan)
% shadedErrorBar(1:size(spkWav,2),mean(spkWav),std(spkWav),'b')
% axis tight
% title(['chan ' num2str(ichan)])
% end
%  
% ichan = 9;
% load([recPath 'spikes\spk_' num2str(ichan)]);
% spkTime(spkTime<1) = [];
% spkTime(spkTime>(size(lfpMat,2)/lfpFs)-1) = [];
%  
% figure(); hold on
% win = -500:500;
% erp = nan(numel(spkTime),numel(win));
% for ispk = 1:numel(spkTime)
%     spkSamp = round(spkTime(ispk)*lfpFs);
%     erp(ispk,:) = lfpMat(31,spkSamp+win);
% end
% sta = mean(erp);
% plot(win/1000,sta,'b')
% hold on
%  
% win = -500:500;
% erp = nan(numel(spkTime),numel(win));
% for ispk = 1:numel(spkTime)
%     spkSamp = round(spkTime(ispk)*lfpFs);
%     erp(ispk,:) = lfpMat(32+31,spkSamp+win);
% end
% sta = mean(erp);
% plot(win/1000,sta,'r')
% hold on
%  
% win = -500:500;
% erp = nan(numel(spkTime),numel(win));
% for ispk = 1:numel(spkTime)
%     spkSamp = round(spkTime(ispk)*lfpFs);
%     erp(ispk,:) = lfpMat(64+4,spkSamp+win);
% end
% sta = mean(erp);
% plot(win/1000,sta,'k')
% hold on
%  
% legend('PPC LFP','VC LFP','LP/Pulvinar LFP')

%% plot coherence couple to LPl phase

lowFreq      = 1; %0.5
highFreq     = 80; %128
numFreqs     = 60; %80
lfpFs = 1000;
foi          = logspace(log10(lowFreq),log10(highFreq),numFreqs);
wavs         = is_makeWavelet(foi,lfpFs);


% compute the complex component of PPC, VC and pulvinar using wavelets
ppc_c = nan(numFreqs,size(lfpMat,2));
vc_c  = nan(numFreqs,size(lfpMat,2));
pul_c = nan(numFreqs,size(lfpMat,2));
selPPC_lfp = lfpMat(1,:);%(0+2,:);
selVC_lfp  = lfpMat(3,:);%(65+14,:);
selPUL_lfp = lfpMat(2,:);%(41+0,:);

ptic %~1.5min
for f = 1:numFreqs
    display(num2str(f))
    ppc_c(f,:) = conv(selPPC_lfp,wavs{f},'same');
    vc_c(f,:)  = conv(selVC_lfp,wavs{f},'same');
    pul_c(f,:) = conv(selPUL_lfp,wavs{f},'same');
end
ptoc


% compute coherence between PPC and VC
Sxy = nanmean(ppc_c.*conj(vc_c),2);
Sxx = nanmean(ppc_c.*conj(ppc_c),2);
Syy = nanmean(vc_c.*conj(vc_c),2);
Cy  = Sxy./(sqrt(Sxx.*Syy));

doPlot = 0;
if doPlot == 1
    figure() 
    plot(foi,abs(Cy),'b','LineWidth',2)
    xlabel('LFP Frequency')
    ylabel('Coherence')
    title('Visual cortex to PPC')

    fois = [0.5 1 2 4 8 16 32 64 128];
    for fi = 1:numel(fois)
        [bi,bb] = sort(abs(foi-fois(fi)));
        tickLoc(fi) = bb(1);
    end
    tickLabel = {'0.5','1','2','4','8','16','32','64','128'};

    figure()
    plot(1:numel(foi),abs(Cy),'b','LineWidth',2)
    set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
    xlabel('LFP Frequency')
    ylabel('Coherence')
    title('Visual cortex to PPC functional connectivity')


    % plot LP/Pulvinar power spectrum
    fig = figure();
    pulPow = nanmean(abs(pul_c).^2,2);
    loglog(foi,pulPow,'LineWidth',2)
    title('LP/Pulvinar power spectrum')
    xlabel('LFP Frequency')
    ylabel('power')
    savefig(fig, 'LPl power spectrum.fig');
end

%% Calculate coherence coupling to selected LPl freqs -- 
% selected frequencies: f=34=LPl peak theta(5Hz), 44=LPl peak alpha (10.23Hz), 49=0153 endogenous alpha(13.5Hz)
% plotCorticalCohCoupleLPl(34, 5);
% plotCorticalCohCoupleLPl(44, 10);
% plotCorticalCohCoupleLPl(49, 13);

% All LPl frequencies
LPl_f = 1:0.5:20;
ksresult = [];
for iLPl_f = 1:numel(LPl_f)
    realf = LPl_f(iLPl_f);
    [bb,bi] = sort(abs(foi-realf)); % find the foi closet to the LPl freq
    cohHist(iLPl_f,:) = plotCorticalCohCoupleLPl(bi(1), realf,pul_c, ppc_c, vc_c,foi,saveRootDir);
    [h,p] = kstest(cohHist(iLPl_f,:));
    ksresult = [ksresult; h p];
end

nbin       = 589;
phaseBins  = 0:2*pi/nbin:2*pi;

save([saveRootDir 'cohHist.mat'], 'cohHist','LPl_f','phaseBins', 'ksresult');

% Plot cohHist for each LPl frequency
screensize = get(groot, 'Screensize' );
fig = figure('Position',[10 50 screensize(3) (screensize(4)-150)/1.2]);
subplot(1,2,1)
imagesc(phaseBins,LPl_f,cohHist);
set(gca,'TickDir','out','YDir','normal','XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})  % to inverse Y axis
ylabel('LPl frequency [Hz]')
xlabel('LPl phase frequency [Hz]')
colormap(jet); cl = colorbar; ylabel(cl,'PPC/V1 gamma coherence (30-80Hz)','FontSize',12); 
caxis([0.1,0.35]);

subplot(1,2,2) % smoothed figure
pcolor(phaseBins,LPl_f,cohHist);
shading interp;
set(gca,'TickDir','out','YDir','normal','XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})  % to inverse Y axis
ylabel('LPl frequency [Hz]')
xlabel('LPl phase frequency [Hz]')
colormap(jet); cl = colorbar; ylabel(cl,'PPC/V1 gamma coherence (30-80Hz)','FontSize',12); 
caxis([0.1,0.35]);
ylim([2,20]);

% subplot(1,2,2) % plot ksresult -- all significant
% plot(LPl_f,ksresult(:,2)); % plot p values for each frequency
% ylabel('KS test p-value')
% xlabel('LPl phase frequency [Hz]')
% view([-90 90])
savefig(fig, [saveRootDir 'figures\cortical gamma_LPl_all_30-60avg.fig']);
saveas(fig, [saveRootDir 'figures\cortical gamma_LPl_all_30-60avg.png']);

%%


function cohHist = plotCorticalCohCoupleLPl(f, realf, pul_c, ppc_c, vc_c,foi,saveRootDir)

ang = mod(angle(pul_c(f,:)),2*pi); % (1xtime) get phase of LPl at a selected frequency
nbin       = 589;
phaseBins  = 0:2*pi/nbin:2*pi;
winSz = pi/8; % window size
winSt = 0.01; % window step
winCen = winSz/2:winSt:2*pi-winSz/2;
ptic
for ibin = 1:numel(winCen)
    display([num2str(ibin) '/' num2str(numel(winCen))])
    pb1 = winCen(ibin)-winSz/2;
    pb2 = winCen(ibin)+winSz/2;
    samps = (ang>=pb1&ang<=pb2); % get time columns where LPl phase is within phase window
    Sxy = nanmean(ppc_c(:,samps).*conj(vc_c(:,samps)),2); % get coherence of these time columns
    Sxx = nanmean(ppc_c(:,samps).*conj(ppc_c(:,samps)),2);
    Syy = nanmean(vc_c(:,samps).*conj(vc_c(:,samps)),2);
    Cy  = Sxy./(sqrt(Sxx.*Syy));
    phase_cy(ibin,:) = Cy;
end
ptoc
mat = flipud(rot90(abs(phase_cy)));
lowFreq      = 0.5;
highFreq     = 128;
numFreqs     = 1000;
intF         = linspace(lowFreq,highFreq,numFreqs);
X = repmat(winCen,numel(foi),1);
Xv = repmat(winCen,numel(intF),1);
Y = repmat(foi',1,numel(winCen));
Yv = repmat(intF',1,numel(winCen));
intMat = interp2(X,Y,mat,Xv,Yv);

f_low = 30; 
f_high = 60;
f_low_ind = floor((f_low-lowFreq)/(highFreq-lowFreq)*(numFreqs-1)+1);
f_high_ind = ceil((f_high-lowFreq)/(highFreq-lowFreq)*(numFreqs-1)+1);
cohHist = mean(intMat(f_low_ind:f_high_ind,:),1);


doPlot = 0;
if doPlot == 1
    % plot cortical coherence coupling with LPl phase
    screensize = get(groot, 'Screensize' );
    fig = figure('Position',[10 50 screensize(3)-100 (screensize(4)-150)/2.5]);

    subplot(1,3,1)
    imagesc(phaseBins,intF,intMat);
    ylim([1 80]);
    set(gca,'TickDir','out','YDir','normal','XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
    %colormap(awesomeMap);
    colormap(jet)
    %caxis([0 0.06])
    y = colorbar;
    ylabel(y,'Visual cortex to PPC coherence')
    xlabel('LP/Pulvinar theta phase')
    ylabel('LFP frequency (Hz)')

    subplot(1,3,2)
    imagesc(phaseBins,intF,intMat);
    ylim([30 60])
    set(gca,'TickDir','out','YDir','normal','XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
    %colormap(awesomeMap);
    colormap(jet)
    caxis([0 0.5]);
    y = colorbar;
    ylabel(y,'Visual cortex to PPC coherence')
    xlabel('LP/Pulvinar theta phase')
    ylabel('LFP frequency (Hz)')


    % plot histogram of cortical coherence of a frequency band over phase
    subplot(1,3,3)
    bar(phaseBins, mean(intMat(f_low_ind:f_high_ind,:),1));
    set(gca,'TickDir','out','YDir','normal','XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
    ylabel('Average PPC-VC coherence')
    xlabel('LP/Pulvinar theta phase')
    savefig(fig, [saveRootDir 'figures\cortical gamma_LPl ' num2str(realf) 'Hz_30-60avg.fig']);
    saveas(fig, [saveRootDir 'figures\cortical gamma_LPl ' num2str(realf) 'Hz_30-60avg.png']);
end

% Putting together all LPl frequency from 1-20Hz


end
% log scale
% fois = [0.5 1 2 4 8 16 32 64 128];
% for fi = 1:numel(fois)
%     [bi,bb] = sort(abs(foi-fois(fi)));
%     tickLoc(fi) = bb(1);
% end
% tickLabel = {'0.5','1','2','4','8','16','32','64','128'};
% 
% mat = flipud(rot90(abs(phase_cy)));
% imagesc(0:0.01:2*pi,1:numel(foi),mat)
% set(gca,'YDir','normal','YTick',tickLoc,'YTickLabel',tickLabel,'XTick',[0 pi/2 pi 3*pi/2 2*pi])
% ylim([55 70])
% colormap(awesomeMap)
% caxis([0 0.15])

