%% After getting funcCon structure, load it from different sessions and do stats

ephysDir     = ['J:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\Analyzed\'];     
saveRootPath = ['J:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\5CSRTT\TimeWindowANOVA\0153_AttentionTask3\'];
%load([ephysDir])


%% define factors: time window and freq band
% time windows of interest
timeWin(1,:) = [-4,-2];
timeWin(2,:) = [-2,0];
timeWin(3,:) = [0,2];
t = size(timeWin,1);

% frequency bands
freqBand(1,:) = [4,8];  % theta
freqBand(2,:) = [8,15]; % alpha
freqBand(3,:) = [15,25];% beta
freqBand(4,:) = [40,80];% gamma

PPC_mean = [];
LPl_mean = [];
VC_mean =[];
PPC_LPl_plv = [];
PPC_LPl_coh = [];
PPC_VC_plv = [];
PPC_VC_coh = [];
LPl_VC_plv = [];
LPl_VC_coh = [];

% % normalize power across all freq bands and entire duration 
% xspecNorm = funcCon.xspec ./ sum(sum(funcCon.xspec(:,foi>freqBand(1,1)&foi<=freqBand(4,2))))*100; 
%1-2
%load('E:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\Analyzed\0153_AttentionTask3_11_20171020_171020_142407\funcConPCA_region1-2.mat')
%load('E:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\Analyzed\0153_AttentionTask3_12_20171024_171024_164233\funcConPCA_region1-2.mat')
%load('E:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\Analyzed\0153_AttentionTask3_13_20171025_171025_125607\funcConPCA_region1-2.mat')
load('E:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\Analyzed\0153_AttentionTask3_14_20171027_171027_171016\funcConPCA_region1-2.mat')

[x_meanPower,x_maxFreq] = TimeWindowMean(freqBand,timeWin,funcCon.tvec,funcCon.xspec);
[y_meanPower,y_maxFreq] = TimeWindowMean(freqBand,timeWin,funcCon.tvec,funcCon.yspec);
[plv_meanPower,plv_maxFreq] = TimeWindowMean(freqBand,timeWin,funcCon.tvec,funcCon.plv);
[coh_meanPower,coh_maxFreq] = TimeWindowMean(freqBand,timeWin,funcCon.tvec,funcCon.coherence);

PPC_mean = [PPC_mean; x_meanPower x_maxFreq];
LPl_mean = [LPl_mean; y_meanPower y_maxFreq];
PPC_LPl_plv = [PPC_LPl_plv; plv_meanPower,plv_maxFreq];
PPC_LPl_coh = [PPC_LPl_coh; coh_meanPower,coh_maxFreq];

%1-3
%load('E:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\Analyzed\0153_AttentionTask3_11_20171020_171020_142407\funcConPCA_region1-3.mat')
%load('E:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\Analyzed\0153_AttentionTask3_12_20171024_171024_164233\funcConPCA_region1-3.mat')
%load('E:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\Analyzed\0153_AttentionTask3_13_20171025_171025_125607\funcConPCA_region1-3.mat')
load('E:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\Analyzed\0153_AttentionTask3_14_20171027_171027_171016\funcConPCA_region1-3.mat')

[y_meanPower,y_maxFreq] = TimeWindowMean(freqBand,timeWin,funcCon.tvec,funcCon.yspec);
[plv_meanPower,plv_maxFreq] = TimeWindowMean(freqBand,timeWin,funcCon.tvec,funcCon.plv);
[coh_meanPower,coh_maxFreq] = TimeWindowMean(freqBand,timeWin,funcCon.tvec,funcCon.coherence);
VC_mean = [VC_mean; y_meanPower y_maxFreq];
PPC_VC_plv = [PPC_VC_plv; plv_meanPower,plv_maxFreq];
PPC_VC_coh = [PPC_VC_coh; coh_meanPower,coh_maxFreq];

%2-3
%load('E:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\Analyzed\0153_AttentionTask3_11_20171020_171020_142407\funcConPCA_region2-3.mat')
%load('E:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\Analyzed\0153_AttentionTask3_12_20171024_171024_164233\funcConPCA_region2-3.mat')
%load('E:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\Analyzed\0153_AttentionTask3_13_20171025_171025_125607\funcConPCA_region2-3.mat')
load('E:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\Analyzed\0153_AttentionTask3_14_20171027_171027_171016\funcConPCA_region2-3.mat')

[plv_meanPower,plv_maxFreq] = TimeWindowMean(freqBand,timeWin,funcCon.tvec,funcCon.plv);
[coh_meanPower,coh_maxFreq] = TimeWindowMean(freqBand,timeWin,funcCon.tvec,funcCon.coherence);
LPl_VC_plv = [LPl_VC_plv; plv_meanPower,plv_maxFreq];
LPl_VC_coh = [LPl_VC_coh; coh_meanPower,coh_maxFreq];

save('0153_14.mat','PPC_mean','LPl_mean','VC_mean','PPC_LPl_plv','PPC_LPl_coh','PPC_VC_plv','PPC_VC_coh','LPl_VC_plv','LPl_VC_coh', '-v7.3');

Allmean = [PPC_mean(:,1) PPC_mean(:,4) LPl_mean(:,1) LPl_mean(:,4) VC_mean(:,1) VC_mean(:,4)...
    PPC_LPl_plv(:,1) PPC_LPl_plv(:,4) PPC_LPl_coh(:,1) PPC_LPl_coh(:,4)...
    PPC_VC_plv(:,1) PPC_VC_plv(:,4) PPC_VC_coh(:,1) PPC_VC_coh(:,4)...
    LPl_VC_plv(:,1) LPl_VC_plv(:,4) LPl_VC_coh(:,1) LPl_VC_coh(:,4);...
    PPC_mean(:,2) PPC_mean(:,5) LPl_mean(:,2) LPl_mean(:,5) VC_mean(:,2) VC_mean(:,5)...
    PPC_LPl_plv(:,2) PPC_LPl_plv(:,5) PPC_LPl_coh(:,2) PPC_LPl_coh(:,5)...
    PPC_VC_plv(:,2) PPC_VC_plv(:,5) PPC_VC_coh(:,2) PPC_VC_coh(:,5)...
    LPl_VC_plv(:,2) LPl_VC_plv(:,5) LPl_VC_coh(:,2) LPl_VC_coh(:,5);...
    PPC_mean(:,3) PPC_mean(:,6) LPl_mean(:,3) LPl_mean(:,6) VC_mean(:,3) VC_mean(:,6)...
    PPC_LPl_plv(:,3) PPC_LPl_plv(:,6) PPC_LPl_coh(:,3) PPC_LPl_coh(:,6)...
    PPC_VC_plv(:,3) PPC_VC_plv(:,6) PPC_VC_coh(:,3) PPC_VC_coh(:,6)...
    LPl_VC_plv(:,3) LPl_VC_plv(:,6) LPl_VC_coh(:,3) LPl_VC_coh(:,6)];

save('Allmean.mat','Allmean');

% region_theta_pow = [PPC_mean(1:4:end,1:t) LPl_mean(1:4:end,1:t) VC_mean(1:4:end,1:t)]; %3 time region
% region_alpha_pow = [PPC_mean(2:4:end,1:t) LPl_mean(2:4:end,1:t) VC_mean(2:4:end,1:t)];
% region_beta_pow  = [PPC_mean(3:4:end,1:t) LPl_mean(3:4:end,1:t) VC_mean(3:4:end,1:t)];
% region_gamma_pow = [PPC_mean(4:4:end,1:t) LPl_mean(4:4:end,1:t) VC_mean(4:4:end,1:t)];
% 
% region_theta_f = [PPC_mean(1:4:end,t+1:end) LPl_mean(1:4:end,t+1:end) VC_mean(1:4:end,t+1:end)]; %3 time region
% region_alpha_f = [PPC_mean(2:4:end,t+1:end) LPl_mean(2:4:end,t+1:end) VC_mean(2:4:end,t+1:end)];
% region_beta_f  = [PPC_mean(3:4:end,t+1:end) LPl_mean(3:4:end,t+1:end) VC_mean(3:4:end,t+1:end)];
% region_gamma_f = [PPC_mean(4:4:end,t+1:end) LPl_mean(4:4:end,t+1:end) VC_mean(4:4:end,t+1:end)];
% 
% plv_theta = [PPC_LPl_plv(1:4:end,:); PPC_VC_plv(1:4:end,:); LPl_VC_plv(1:4:end,:)];
% plv_alpha = [PPC_LPl_plv(2:4:end,:); PPC_VC_plv(2:4:end,:); LPl_VC_plv(2:4:end,:)];
% plv_beta  = [PPC_LPl_plv(3:4:end,:); PPC_VC_plv(3:4:end,:); LPl_VC_plv(3:4:end,:)];
% plv_gamma = [PPC_LPl_plv(4:4:end,:); PPC_VC_plv(4:4:end,:); LPl_VC_plv(4:4:end,:)];




doPlot = 1;
if doPlot == 1
    screensize = get( groot, 'Screensize' );
    fig = figure('Position',[10 50 screensize(3)-100 screensize(4)-150]);
    subplot(2,1,1)
    bar(pow2db(x_meanPower));
    ylabel('Power [dB] signal X','FontSize',12);xlabel('Frequency band');
    
    subplot(2,1,2)
    bar(pow2db(y_meanPower));
    ylabel('Power [dB] signal Y','FontSize',12);xlabel('Frequency band');
    savefig(fig,['figures\level3_11_PPC&LPl_meanPower.fig'],'compact');
    
    fig = figure('Position',[10 50 screensize(3)-100 screensize(4)-150]);
    subplot(2,1,1)
    bar(x_maxFreq);
    ylabel('Peak power frequency [Hz]','FontSize',12);xlabel('Frequency band');
    set(gca,'XTickLabel',{'Theta (4-8Hz)' 'Alpha (8-15Hz)','Beta (14-25Hz)','Gamma (25-50Hz)'});

    
    subplot(2,1,2)
    bar(y_maxFreq);
    ylabel('Peak power frequency [Hz]','FontSize',12);xlabel('Frequency band');
    set(gca,'XTickLabel',{'Theta (4-8Hz)' 'Alpha (8-14Hz)','Beta (14-25Hz)','Gamma (25-50Hz)'});
    savefig(fig,[saveRootPath 'figures\level3_11_PPC&LPl_maxFreq.fig'],'compact');
    
end



