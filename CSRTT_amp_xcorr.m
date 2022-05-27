clear

cluster = 0;
tic

%% load data and functions

if cluster == 0
    addpath('E:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys'); % all the functions needed
    addpath('E:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\mvgc_v1.0') % for Granger Causality
    ephysDir     = ['E:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\toAnalyze\0153_AttentionTask3_13\'];     
    saveRootPath = ['E:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\0153\Analyzed\0153_AttentionTask3_13\amp_xcorr\']; 
    mkdir(saveRootPath);  
    %addpath(genpath(ephysDir)); % add folder and subfolders to path
    cd(saveRootPath); % change directory to save data

    % load behavioral data
    load('E:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\5CSRTT\Behav data\0153\0153_Level3_13_20171025.mat'); %20,24,25,27
    % load lfpData
    load([ephysDir 'lfp\lfpPCA.mat']);
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
fs = lfpFs;
twin = [-5 8];
onInitInd = 1;
stimTouchInd = 2;

rawFs = 30000;
trialOnset = find(diff(triggerData(onInitInd,:))==1)./rawFs;
trialInit = find(diff(triggerData(onInitInd,:))==-1)./rawFs;
stimOnset = find(diff(triggerData(stimTouchInd,:))==1)./rawFs;
touch = find(diff(triggerData(stimTouchInd,:))==-1)./rawFs;

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

%%
% for 0153_level3: PPC, LPl, VC = 1,2,3 row of lfpPCA
f_theta = [3,6];
low_freq = f_theta(1);
high_freq = f_theta(2);

segLength = 7; % in sec, needs to be at least 3 for a filter with order=1000 (in amp_xcorr_3region.m) 
newFs     = 400;  
maxLag = 4*round(newFs/10); %round(fs/10)
numLags = 2*maxLag +1;

xsig = lfpPCA(1,:);
ysig = lfpPCA(2,:);
zsig = lfpPCA(3,:);

% % reject noise in recording (high amplitude noise can destroy coherence estimates in particular)
% xsig = sub_rejectNoise(xser,fs,2,1); 
% ysig = sub_rejectNoise(yser,fs,2,1); 
% zsig = sub_rejectNoise(zser,fs,2,1); 

% define low pass filter at 100Hz
nyqFreq = fs/2;
[b,a]   = butter(2,100/nyqFreq,'low');
xfilt   = filtfilt(b,a,xsig);
yfilt   = filtfilt(b,a,ysig);
zfilt   = filtfilt(b,a,zsig);

% resample data to sample rate of newFs (defined above)
idat = resample(xfilt,newFs,fs);
jdat = resample(yfilt,newFs,fs);
zdat = resample(zfilt,newFs,fs);

% to get time resolution of Xcorr, divide window of interests into 3 sec
% sliding windows, and calculate for each sliding window 
stepSize  = 0.1; % sliding window increments
stepCen   = twin(1):stepSize:twin(2); % sliding window centers for each calculation
recLength = numel(xsig)/newFs; % whole recording length in sec
halfWinSamp = (segLength*newFs)/2; 

% % power spectrum of the whole recording
% figure()
% ipxx = pwelch(xfilt,newFs);
% jpxx = pwelch(yfilt,newFs);
% zpxx = pwelch(zfilt,newFs);
% plot(10*log10(ipxx));hold on;
% plot(10*log10(jpxx))
% plot(10*log10(zpxx))
% legend('PPC','LPl','V1')

c = 0;% counting variable
for iev = 1:numel(event) % loop through each trial
       
    % skip if we have window range issues
    if event(iev) < abs(twin(1))+segLength; continue; end
    if event(iev) > recLength-twin(2)-segLength; continue; end
    c = c + 1; % count for trial 
    
    for istep = 1:numel(stepCen) % loop through each sliding time window centers
        
        clear imat jmat zmat % due to parfor

        samp      = round((stepCen(istep)+event(iev))*newFs);

        % Cut out data around events    
        imat = idat(samp-halfWinSamp:samp+halfWinSamp);
        jmat = jdat(samp-halfWinSamp:samp+halfWinSamp);
        zmat = zdat(samp-halfWinSamp:samp+halfWinSamp);

        [lags, crosscorr_12, max_crosscorr_lag_12,crosscorr_23, max_crosscorr_lag_23,crosscorr_13, max_crosscorr_lag_13]...
        =amp_xcorr_3region(imat,jmat,zmat,newFs,low_freq,high_freq,maxLag);

        Xcorr.X2Ylag(c,istep) = max_crosscorr_lag_12;
        Xcorr.Y2Zlag(c,istep) = max_crosscorr_lag_23;
        Xcorr.X2Zlag(c,istep) = max_crosscorr_lag_13;        
        Xcorr.X2Ymax(c,istep) = max(crosscorr_12);
        Xcorr.Y2Zmax(c,istep) = max(crosscorr_23);
        Xcorr.X2Zmax(c,istep) = max(crosscorr_13);
                
        Xcorr.X2Ycorr_vec(c,istep,:) = crosscorr_12;
        Xcorr.Y2Zcorr_vec(c,istep,:) = crosscorr_23;
        Xcorr.X2Zcorr_vec(c,istep,:) = crosscorr_13;          
       
        
    end
end

Xcorr.lags = lags;
Xcorr.tvec = stepCen;

% replace all outliers and 0 by nan
Xcorr.X2Ylag_nan = Xcorr.X2Ylag;
Xcorr.X2Ylag_nan(Xcorr.X2Ylag>=200 | Xcorr.X2Ylag<=-200 | Xcorr.X2Ylag== 0) = nan;
Xcorr.Y2Zlag_nan = Xcorr.Y2Zlag;
Xcorr.Y2Zlag_nan(Xcorr.Y2Zlag>=200 | Xcorr.Y2Zlag<=-200 | Xcorr.Y2Zlag== 0) = nan;
Xcorr.X2Zlag_nan = Xcorr.X2Zlag;
Xcorr.X2Zlag_nan(Xcorr.X2Zlag>=200 | Xcorr.X2Zlag<=-200 | Xcorr.X2Zlag== 0) = nan;
save([saveRootPath 'amp_Xcorr_-400~400ms_7sec.mat'],'Xcorr','-v7.3');

toc;

doPlot = 1;
if doPlot == 1
    % plot one trial of xcorr over time

    for itrial = 1:size(Xcorr.X2Ylag,1)
        screensize = get( groot, 'Screensize');
        fig1 = figure('Position',[10 50 (screensize(3)-100)/2 screensize(4)-150]);
        subplot(311)
        imagesc(Xcorr.tvec, Xcorr.lags, squeeze(Xcorr.X2Ycorr_vec(itrial,:,:))'); hold on;
        hline(0,'k-');
        scatter(Xcorr.tvec, Xcorr.X2Ylag(itrial,:),20,'w','filled');
        xlabel('Time to touch [s]'); ylabel('Lags [ms]'); title('PPC-LPl lag')
        set(gca,'YDir','normal','TickDir','out');colormap jet
        cl = colorbar('eastoutside'); ylabel(cl,'X-correlation','FontSize',12)

        subplot(312)
        imagesc(Xcorr.tvec, Xcorr.lags, squeeze(Xcorr.Y2Zcorr_vec(itrial,:,:))'); hold on;
        hline(0,'k-');
        scatter(Xcorr.tvec, Xcorr.Y2Zlag(itrial,:),20,'w','filled');
        xlabel('Time to touch [s]'); ylabel('Lags [ms]'); title('LPl-V1 lag')
        set(gca,'YDir','normal','TickDir','out');colormap jet
        cl = colorbar('eastoutside'); ylabel(cl,'X-correlation','FontSize',12)

        subplot(313)
        imagesc(Xcorr.tvec, Xcorr.lags, squeeze(Xcorr.X2Zcorr_vec(itrial,:,:))'); hold on;
        hline(0,'k-');
        scatter(Xcorr.tvec, Xcorr.X2Zlag(itrial,:),20,'w','filled');
        xlabel('Time to touch [s]'); ylabel('Lags [ms]'); title('PPC-V1 lag')
        set(gca,'YDir','normal','TickDir','out');colormap jet
        cl = colorbar('eastoutside'); ylabel(cl,'X-correlation','FontSize',12)

        if ~exist([saveRootPath 'figures'],'dir'); mkdir([saveRootPath 'figures']); end
        saveas(fig1,[saveRootPath, 'figures/trial_' num2str(itrial) '.png']);
        
        
    end
    
    % box plot of lags across trials for each time window
    for iWindow = 1:size(Xcorr.X2Ylag,2)
        boxPlot = 0;
        if boxPlot == 1
            fig2 = figure('Position',[10 50 (screensize(3)-100)/2 screensize(4)-150]);
            all_lag(:,1) = Xcorr.X2Ylag_nan(:,iWindow);
            all_lag(:,2) = Xcorr.X2Zlag_nan(:,iWindow);
            all_lag(:,3) = Xcorr.Y2Zlag_nan(:,iWindow);
            boxplot(all_lag,'Labels',{'PPC-LPl','LPl-V1', 'PPC-V1'}, 'notch', 'on');
            title(sprintf('PPC PUL=%1.2f, PUL VC=%1.2f, PPC VC=%1.2f',...
                nanmedian(all_lag(:,1)), nanmedian(all_lag(:,2)), ...
                nanmedian(all_lag(:,3))))
            %p = anova1(all_lag);
            saveas(fig2,[saveRootPath, 'figures/box_timeWindow_' num2str(iWindow) '.png']);
        end
        
%         % use histfit can automatically generate fitting curve, but the x
%         % range is not consistent across region and trials, can't use it for heatmap
%         hstCnt = 10;
%         maxLag = 200;
%         hiswidth = 2*maxLag/hstCnt;
%         hstLmt = [-maxLag maxLag];
%         hstCtrs = [-maxLag+hiswidth/2:hiswidth:maxLag-hiswidth/2]; % centers
%         
%         %hHist(1).YData is the histogram, hHist(2) is the density curve
%         fig3 = figure('Position',[10 50 (screensize(3)-100)/2 screensize(4)-150]);
%         subplot(311)
%         hHist1 = histfit(Xcorr.X2Ylag_nan(:,iWindow), hstCnt, 'kernel');
%             xlim(hstLmt);
%             allHist.hist(1,iWindow,:) = hHist1(1).YData;
%             allHist.curve(1,iWindow,:) = hHist1(2).YData;
%             allHist.hist_x = hHist1(1).XData;
%             allHist.curve_x = hHist1(2).XData;
%             x_axis = get(gca,'x');
%         subplot(312)
%         hHist2 = histfit(Xcorr.Y2Zlag_nan(:,iWindow), hstCnt, 'kernel');
%             xlim(hstLmt);
%             allHist.hist(2,iWindow,:) = hHist2(1).YData;
%             allHist.curve(2,iWindow,:) = hHist2(2).YData;            
%         subplot(313)
%         hHist3 = histfit(Xcorr.X2Zlag_nan(:,iWindow), hstCnt, 'kernel');
%             xlim(hstLmt);
%             allHist.hist(3,iWindow,:) = hHist3(1).YData;
%             allHist.curve(3,iWindow,:) = hHist3(2).YData;
%             allHist.hist_x = hHist3(1).XData;
%             allHist.curve_x = hHist3(2).XData;
%         saveas(fig3,[saveRootPath, 'figures/hist_timeWindow_' num2str(iWindow) '.png']);
        
        % use norfit and normpdf to generate normal curve based on data
        lag_interval = 2; % in ms
        lag_x = [-200:lag_interval:200];
        
        temp1 = Xcorr.X2Ylag_nan(:,iWindow);
        temp1 = temp1(~isnan(temp1)); % exclude all nan        
        [muhat,sigmahat] = normfit(temp1);
        norm(1,iWindow,:) = normpdf(lag_x, muhat,sigmahat);

        temp1 = Xcorr.Y2Zlag_nan(:,iWindow);
        temp1 = temp1(~isnan(temp1)); % exclude all nan        
        [muhat,sigmahat] = normfit(temp1);
        norm(2,iWindow,:) = normpdf(lag_x, muhat,sigmahat);

        temp1 = Xcorr.X2Zlag_nan(:,iWindow);
        temp1 = temp1(~isnan(temp1)); % exclude all nan        
        [muhat,sigmahat] = normfit(temp1);
        norm(3,iWindow,:) = normpdf(lag_x, muhat,sigmahat);
        %fig4 = figure('Position',[10 50 (screensize(3)-100)/2 screensize(4)-150]);
        %plot(lag_x, norm(1,iWindow,:), 'linewidth', 3); 
        %plot([muhat muhat], [0, max(norm)],  'linewidth', 2)
            
    end
    
    save([saveRootPath 'lag_norm_-200~200.mat'],'norm','-v7.3');
    
    % plot histogram of all lags
    hstCnt = 100;
    hstLmt = [-200 200];
    figure()
    subplot(311);
    hHist_12 = histfit(Xcorr.X2Ylag_nan(:), hstCnt,'normal');        
    title('PPC-LPl lag, theta'); xlim(hstLmt);
    xlabel('Lags [ms]');
    subplot(312);
    hHist_23 = histfit(Xcorr.Y2Zlag_nan(:), hstCnt,'normal');
    title('LPl-V1 lag, theta'); xlim(hstLmt)
    subplot(313);
    hHist_13 = histfit(Xcorr.X2Zlag_nan(:), hstCnt,'normal');
    title('PPC-V1 lag, theta'); xlim(hstLmt)
    xlabel('Lags [ms]');
    
    % lag normed distribution heatmap over time
    colorLim = [1e-3, 5e-3];
    screensize = get( groot, 'Screensize');
    fig1 = figure('Position',[10 50 (screensize(3)-100)/2 screensize(4)-150]);
    subplot(311)
        imagesc(Xcorr.tvec, lag_x, squeeze(norm(1,:,:))'); hold on;
        hline(0,'k-');
        scatter(Xcorr.tvec, nanmedian(Xcorr.X2Ylag_nan,1),20,'w','filled');
        xlabel('Time to touch [s]'); ylabel('Lags [ms]'); title('PPC-LPl lag')
        set(gca,'YDir','normal','TickDir','out');colormap jet
        cl = colorbar('eastoutside'); ylabel(cl,'Lag normed-distribution','FontSize',12)
        caxis(colorLim);
        
    subplot(312)
        imagesc(Xcorr.tvec, lag_x, squeeze(norm(2,:,:))'); hold on;
        hline(0,'k-');
        scatter(Xcorr.tvec, nanmedian(Xcorr.Y2Zlag_nan,1),20,'w','filled');
        xlabel('Time to touch [s]'); ylabel('Lags [ms]'); title('LPl-V1 lag')
        set(gca,'YDir','normal','TickDir','out');colormap jet
        cl = colorbar('eastoutside'); ylabel(cl,'Lag normed-distribution','FontSize',12)
        caxis(colorLim);

    subplot(313)
        imagesc(Xcorr.tvec, lag_x, squeeze(norm(3,:,:))'); hold on;
        hline(0,'k-');
        scatter(Xcorr.tvec, nanmedian(Xcorr.X2Zlag_nan,1),20,'w','filled');
        xlabel('Time to touch [s]'); ylabel('Lags [ms]'); title('PPC-V1 lag')
        set(gca,'YDir','normal','TickDir','out');colormap jet
        cl = colorbar('eastoutside'); ylabel(cl,'Lag normed-distribution','FontSize',12)
        caxis(colorLim);

        savefig(fig1,[saveRootPath, 'lag_density_over_time.fig'],'compact');
        saveas(fig1,[saveRootPath, 'lag_density_over_time.png']);
         
       
%     % plot histogram heatmap over time (not continuous, not using)
%     screensize = get( groot, 'Screensize');
%     fig1 = figure('Position',[10 50 (screensize(3)-100)/2 screensize(4)-150]);
%     subplot(311)
%         imagesc(Xcorr.tvec, allHist.x, squeeze(allHist.curve(1,:,:))'); hold on;
%         hline(0,'k-');
%         scatter(Xcorr.tvec, nanmedian(Xcorr.X2Ylag_nan,1),20,'w','filled');
%         xlabel('Time to touch [s]'); ylabel('Lags [ms]'); title('PPC-LPl lag')
%         set(gca,'YDir','normal','TickDir','out');colormap jet
%         cl = colorbar('eastoutside'); ylabel(cl,'Lag density','FontSize',12);
% 
%     subplot(312)
%         imagesc(Xcorr.tvec, allHist.x, squeeze(allHist.curve(2,:,:))'); hold on;
%         hline(0,'k-');
%         scatter(Xcorr.tvec, nanmedian(Xcorr.Y2Zlag_nan,1),20,'w','filled');
%         xlabel('Time to touch [s]'); ylabel('Lags [ms]'); title('LPl-V1 lag')
%         set(gca,'YDir','normal','TickDir','out');colormap jet
%         cl = colorbar('eastoutside'); ylabel(cl,'X-correlation','FontSize',12)
% 
%     subplot(313)
%         imagesc(Xcorr.tvec, allHist.x, squeeze(allHist.curve(3,:,:))'); hold on;
%         hline(0,'k-');
%         scatter(Xcorr.tvec, nanmedian(Xcorr.X2Zlag_nan,1),20,'w','filled');
%         xlabel('Time to touch [s]'); ylabel('Lags [ms]'); title('PPC-V1 lag')
%         set(gca,'YDir','normal','TickDir','out');colormap jet
%         cl = colorbar('eastoutside'); ylabel(cl,'X-correlation','FontSize',12)
% 
%         if ~exist([saveRootPath 'figures'],'dir'); mkdir([saveRootPath 'figures']); end
%         saveas(fig2,[saveRootPath, 'lag_histogram_over_time.fig', fig]);
%         saveas(fig2,[saveRootPath, 'lag_histogram_over_time.png']);
    
    % plot median    
    figure()
    subplot(311); plot(Xcorr.tvec,median(Xcorr.X2Ylag(1,:,:),1));title('PPC-LPl lag')
    subplot(312); plot(Xcorr.tvec,median(Xcorr.Y2Zlag(1,:,:),1));title('LPl-V1 lag')
    subplot(313); plot(Xcorr.tvec,median(Xcorr.X2Zlag(1,:,:),1));title('PPC-V1 lag')
    xlabel('Time to touch [s]') 
    
    
    % plot lag heatmap over time and histogram for all trials (combine
    % above 2 figures)
    screensize = get( groot, 'Screensize');
    fig = figure('Position',[10 50 (screensize(3)-100)/2 (screensize(4)-150)/2]);
   
    subplot(3,4,[1 2 3]); 
        imagesc(Xcorr.tvec, lag_x, squeeze(norm(1,:,:))'); hold on;
        hline(0,'k-');
        scatter(Xcorr.tvec, nanmedian(Xcorr.X2Ylag_nan,1),20,'w','filled');
        %xlabel('Time to touch [s]'); 
        ylabel('Lags [ms]'); title('PPC-LPl lag')
        set(gca,'YDir','normal','TickDir','out');colormap jet
        cl = colorbar('eastoutside'); %ylabel(cl,'Lag normed-distribution','FontSize',12)
        caxis(colorLim);
        
    subplot(3,4,[5 6 7])
        imagesc(Xcorr.tvec, lag_x, squeeze(norm(2,:,:))'); hold on;
        hline(0,'k-');
        scatter(Xcorr.tvec, nanmedian(Xcorr.Y2Zlag_nan,1),20,'w','filled');
        %xlabel('Time to touch [s]'); 
        ylabel('Lags [ms]'); title('LPl-V1 lag')
        set(gca,'YDir','normal','TickDir','out');colormap jet
        cl = colorbar('eastoutside'); ylabel(cl,'Lag normed-distribution','FontSize',12)
        caxis(colorLim);

    subplot(3,4,[9 10,11])
        imagesc(Xcorr.tvec, lag_x, squeeze(norm(3,:,:))'); hold on;
        hline(0,'k-');
        scatter(Xcorr.tvec, nanmedian(Xcorr.X2Zlag_nan,1),20,'w','filled');
        xlabel('Time to touch [s]'); ylabel('Lags [ms]'); title('PPC-V1 lag')
        set(gca,'YDir','normal','TickDir','out');colormap jet
        cl = colorbar('eastoutside'); %ylabel(cl,'Lag normed-distribution','FontSize',12)
        caxis(colorLim);

    subplot(3,4,4);
        hstCnt = 40;
        hstLmt = [-200 200];
        yLmt = [0,150];
        histfit(Xcorr.X2Ylag_nan(:), hstCnt,'normal');        
        x = median(nanmedian(Xcorr.X2Ylag_nan,1));
        vline(x,'r--');    
        xlim(hstLmt);ylim(yLmt);
        text(x,0,num2str(x),'color','red','HorizontalAlignment','right','VerticalAlignment','middle');
        view([90 90]);set(gca,'XDir','reverse');
        title('Lag distribution');
        
    subplot(3,4,8);
        histfit(Xcorr.Y2Zlag_nan(:), hstCnt,'normal');        
        x = median(nanmedian(Xcorr.Y2Zlag_nan,1));
        vline(x,'r--');    
        xlim(hstLmt);ylim(yLmt);
        text(x,0,num2str(x),'color','red','HorizontalAlignment','right','VerticalAlignment','middle');
        view([90 90]);set(gca,'XDir','reverse');
        
    subplot(3,4,12);
        histfit(Xcorr.X2Zlag_nan(:), hstCnt,'normal');        
        x = median(nanmedian(Xcorr.X2Zlag_nan,1));
        vline(x,'r--');    
        xlim(hstLmt);ylim(yLmt);
        text(x,0,num2str(x),'color','red','HorizontalAlignment','right','VerticalAlignment','middle');
        view([90 90]);set(gca,'XDir','reverse');

    savefig(fig,[saveRootPath, 'lag_density_over_time+histo.fig'],'compact');
    saveas(fig,[saveRootPath, 'lag_density_over_time+histo.png']);
    
    


end   


