animals = {'0114','0116','0124','0125'};
pathDir = 'D:\';

spkWin = 0.5;

count    = 0;
picLat   = [];
picOnPup = [];
vidLat   = [];
vidOnPup = [];
picPup_ppc  = [];
vidPup_ppc  = [];
picPup_pul  = [];
vidPup_pul  = [];
ppcSpkPic = [];
ppcSpkVid = [];
pulSpkPic = [];
pulSpkVid = [];
for ianimal = 1:numel(animals)
    animalCode    = animals{ianimal};
    socialDarkDir = dir([pathDir animalCode '\' animalCode '_imagesVideo*']);
    
    for irec = 1:numel(socialDarkDir);
        count = count + 1;
        an(count) = ianimal;
        recName = socialDarkDir(irec).name;
        recPath = [pathDir animalCode '\' recName '\'];
        display(recName)
        % recName = '0125_imagesVideo_8_160309';
        % recPath = ['E:\0125\' recName '\'];
        load([recPath 'adc_data'])
        load([recPath 'triggerData'])
        lfpFs = is_load([recPath 'lfp\lfpMat'],'lfpFs');

        % cz_detectSaccadesEphys(recPath,animalCode);
        load([recPath 'pupil'])
        load([recPath 'saccades'])
        mu     = nanmean(pupClean);
        sigma  = nanstd(pupClean);
        sigma0 = sigma;
        sigma0(sigma0==0) = 1;
        z      = bsxfun(@minus,pupClean, mu);
        zpup   = bsxfun(@rdivide, z, sigma0);
        
        try
            % read in the log file
            fileName = ['D:\Presentation_LogFiles\' recName(1:end-7) '.log'];
            fileID   = fopen(fileName);
            formatSpec = '%f %s %s %f %f %f %f %f %f %f %s %f';
            LogFile = textscan(fileID,formatSpec,'HeaderLines',5,'Delimiter', '\t');
            fclose(fileID);
        catch; continue; end
        
        % detect images
        cellPic = cellfun(@regexp,LogFile(2),{'Picture'},'UniformOutput',false);
        picInd  = find(cellfun(@isempty,cellPic{:}) == 0);
        % detect videos
        cellVid = cellfun(@regexp,LogFile(2),{'Video'},'UniformOutput',false);
        vidInd  = find(cellfun(@isempty,cellVid{:}) == 0);
        
        % detect stim onsets
        onset  = find(diff(triggerData(1,:))==1);
        offset = find(diff(triggerData(2,:))==1);
        
        picOnset = round((onset(picInd)/Fs)*lfpFs);
        vidOnset = round((onset(vidInd)/Fs)*lfpFs);
        
        win = round([-5 15]);
        
        % determine spikes
        spkFiles = dir([recPath 'spikes\spk*']);
        ppcChans = 1:32;
        if numel(spkFiles)>48; pulChans = 41:56; else pulChans = 33:48; end
        
        
        % Do PPC
        for k = 1:numel(ppcChans)
            spk = is_load([recPath 'spikes\cleanSpk_' num2str(ppcChans(k))],'index');
            % spk = spk/1000;
            [tPSTH,PSTH] = is_PSTHstats(picOnset/lfpFs,spk,win);
            ppc_pic(k,:) = PSTH;
            [tPSTH,PSTH] = is_PSTHstats(vidOnset/lfpFs,spk,win);
            ppc_vid(k,:) = PSTH;
            
            % detect number of spikes for each trial in first 1s
            
            % images
            spks = nan(1,numel(picOnset));
            p    = nan(1,numel(picOnset));
            for iev = 1:numel(picOnset)
                evT       = picOnset(iev)/lfpFs;
                spks(iev) = numel(find(spk>evT&spk<evT+spkWin));
                p(iev)    = zpup(picOnset(iev));
            end
            picPup_ppc = [picPup_ppc p];
            ppcSpkPic = [ppcSpkPic zscore(spks)];
            
            % video
            spks = nan(1,numel(vidOnset));
            p    = nan(1,numel(vidOnset));
            for iev = 1:numel(vidOnset)
                evT       = vidOnset(iev)/lfpFs;
                spks(iev) = numel(find(spk>evT&spk<evT+spkWin));
                p(iev)    = zpup(vidOnset(iev));
            end
            vidPup_ppc = [vidPup_ppc p];
            ppcSpkVid = [ppcSpkVid zscore(spks)];
            
        end
        ppcPSTH_pic(count,:) = mean(ppc_pic);
        ppcPSTH_vid(count,:) = mean(ppc_vid);

        % Do Pulvinar
        for k = 1:numel(pulChans)
            spk = is_load([recPath 'spikes\cleanSpk_' num2str(pulChans(k))],'index');
            % spk = spk/1000;
            [tPSTH,PSTH] = is_PSTHstats(picOnset/lfpFs,spk,win);
            pul_pic(k,:) = PSTH;
            [tPSTH,PSTH] = is_PSTHstats(vidOnset/lfpFs,spk,win);
            pul_vid(k,:) = PSTH;    
            
            % detect number of spikes for each trial in first 1s
            
            % images
            spks = nan(1,numel(picOnset));
            p    = nan(1,numel(picOnset));
            for iev = 1:numel(picOnset)
                evT       = picOnset(iev)/lfpFs;
                spks(iev) = numel(find(spk>evT&spk<evT+spkWin));
                p(iev)    = zpup(picOnset(iev));
            end
            picPup_pul = [picPup_pul p];
            pulSpkPic = [pulSpkPic zscore(spks)];
            
            % video
            spks = nan(1,numel(vidOnset));
            p    = nan(1,numel(vidOnset));
            for iev = 1:numel(vidOnset)
                evT       = vidOnset(iev)/lfpFs;
                spks(iev) = numel(find(spk>evT&spk<evT+spkWin));
                p(iev)    = zpup(vidOnset(iev));
            end
            vidPup_pul = [vidPup_pul p];
            pulSpkVid = [pulSpkVid zscore(spks)];                        
        end
        pulPSTH_pic(count,:) = mean(pul_pic);
        pulPSTH_vid(count,:) = mean(pul_vid);
        
    end
end

RdBu=cbrewer('div', 'RdBu', 10, 'linear');

% plot PSTH's
figure
mn_ppc = mean(ppcPSTH_vid);
er_ppc = std(ppcPSTH_vid)/size(ppcPSTH_vid,1);
p_ppc = vertcat(mn_ppc,mn_ppc+er_ppc,mn_ppc-er_ppc)';
mn_pul = mean(pulPSTH_vid);
er_pul = std(pulPSTH_vid)/size(pulPSTH_vid,1);
p_pul = vertcat(mn_pul,mn_pul+er_pul,mn_pul-er_pul)';

[hAx,hLine1,hLine2] = plotyy(tPSTH,p_ppc,tPSTH,p_pul); hold on
set(hLine1(1),'Color',RdBu(1,:))
set(hLine1(2),'Color',RdBu(4,:),'LineStyle','--')
set(hLine1(3),'Color',RdBu(4,:),'LineStyle','--')
set(hLine2(1),'Color',RdBu(10,:))
set(hLine2(2),'Color',RdBu(6,:),'LineStyle','--')
set(hLine2(3),'Color',RdBu(6,:),'LineStyle','--')

set(hAx(1),'YColor',RdBu(2,:),'Ylim',[12 22],'TickDir','out','YTick',[12 14 16 18 20 22])
set(hAx(2),'YColor',RdBu(10,:),'Ylim',[7.5 16],'TickDir','out','YTick',[8 10 12 14 16])
set(get(hAx(1), 'Ylabel'), 'String', 'PPC Firing Rate (Hz)');
set(get(hAx(2), 'Ylabel'), 'String', 'LP/Pulvinar Firing Rate (Hz)');
title('Response to video')


% plot PSTH's
figure
mn_ppc = mean(ppcPSTH_pic);
er_ppc = std(ppcPSTH_pic)/size(ppcPSTH_pic,1);
p_ppc = vertcat(mn_ppc,mn_ppc+er_ppc,mn_ppc-er_ppc)';
mn_pul = mean(pulPSTH_pic);
er_pul = std(pulPSTH_pic)/size(pulPSTH_pic,1);
p_pul = vertcat(mn_pul,mn_pul+er_pul,mn_pul-er_pul)';

[hAx,hLine1,hLine2] = plotyy(tPSTH,p_ppc,tPSTH,p_pul); hold on
set(hLine1(1),'Color',RdBu(1,:))
set(hLine1(2),'Color',RdBu(4,:),'LineStyle','--')
set(hLine1(3),'Color',RdBu(4,:),'LineStyle','--')
set(hLine2(1),'Color',RdBu(10,:))
set(hLine2(2),'Color',RdBu(6,:),'LineStyle','--')
set(hLine2(3),'Color',RdBu(6,:),'LineStyle','--')

set(hAx(1),'YColor',RdBu(2,:),'Ylim',[12 22],'TickDir','out','YTick',[12 14 16 18 20 22])
set(hAx(2),'YColor',RdBu(10,:),'Ylim',[7.5 16],'TickDir','out','YTick',[8 10 12 14 16])
set(get(hAx(1), 'Ylabel'), 'String', 'PPC Firing Rate (Hz)');
set(get(hAx(2), 'Ylabel'), 'String', 'LP/Pulvinar Firing Rate (Hz)');
title('Response to static image')


% do some stats
grabTime = (tPSTH>0 & tPSTH<10);
ppc_pic_av = mean(ppcPSTH_pic(:,grabTime),2);
ppc_vid_av = mean(ppcPSTH_vid(:,grabTime),2);
pul_pic_av = mean(pulPSTH_pic(:,grabTime),2);
pul_vid_av = mean(pulPSTH_vid(:,grabTime),2);

% PPC
mn = mean(ppc_pic_av);
er = std(ppc_pic_av)/sqrt(numel(ppc_pic_av));
display(['PPC pic mean = ' num2str(mn) ' err = ' num2str(er)])
mn = mean(ppc_vid_av);
er = std(ppc_vid_av)/sqrt(numel(ppc_vid_av));
display(['PPC vid mean = ' num2str(mn) ' err = ' num2str(er)])
[h,p] = ttest(ppc_pic_av,ppc_vid_av);
p     = signtest(ppc_pic_av,ppc_vid_av);
display(['PPC stats p value = ' num2str(p)])

% LP/Pulvinar
mn = mean(pul_pic_av);
er = std(pul_pic_av)/sqrt(numel(pul_pic_av));
display(['LP/Pulvinar pic mean = ' num2str(mn) ' err = ' num2str(er)])
mn = mean(pul_vid_av);
er = std(pul_vid_av)/sqrt(numel(pul_vid_av));
display(['LP/Pulvinar vid mean = ' num2str(mn) ' err = ' num2str(er)])
[h,p] = ttest(pul_pic_av,pul_vid_av);
p     = signtest(pul_pic_av,pul_vid_av);
display(['LP/Pulvinar stats p value = ' num2str(p)])



plot(tPSTH,mn+er,'r--')
plot(tPSTH,mn-er,'r--')
xlabel('Time (s)')
ylabel('Spikes/s')
title('PPC')

yyaxis right
plot(tPSTH,mn,'b'); hold on
plot(tPSTH,mn+er,'b--')
plot(tPSTH,mn-er,'b--')

mn = mean(ppcPSTH_pic);
er = std(ppcPSTH_pic)/size(ppcPSTH_pic,1);
plot(tPSTH,mn,'b'); hold on
plot(tPSTH,mn+er,'b--')
plot(tPSTH,mn-er,'b--')




% plot PSTH's
figure
subplot(1,2,1)
mn = mean(ppcPSTH_pic);
er = std(ppcPSTH_pic)/size(ppcPSTH_pic,1);
plot(tPSTH,mn,'b'); hold on
plot(tPSTH,mn+er,'b--')
plot(tPSTH,mn-er,'b--')

mn = mean(ppcPSTH_vid);
er = std(ppcPSTH_vid)/size(ppcPSTH_vid,1);
plot(tPSTH,mn,'r'); hold on
plot(tPSTH,mn+er,'r--')
plot(tPSTH,mn-er,'r--')
xlabel('Time (s)')
ylabel('Spikes/s')
title('PPC')

subplot(1,2,2)
mn = mean(pulPSTH_pic);
er = std(pulPSTH_pic)/size(pulPSTH_pic,1);
plot(tPSTH,mn,'b'); hold on
plot(tPSTH,mn+er,'b--')
plot(tPSTH,mn-er,'b--')

mn = mean(pulPSTH_vid);
er = std(pulPSTH_vid)/size(pulPSTH_vid,1);
plot(tPSTH,mn,'r'); hold on
plot(tPSTH,mn+er,'r--')
plot(tPSTH,mn-er,'r--')
xlabel('Time (s)')
ylabel('Spikes/s')
title('LP/Pulvinar')






% scatter onset spikes vs pupil diameter
a = picPup_ppc;
b = ppcSpkPic;
[t,p] = corrcoef(a,b,'rows','complete')
scatter(a,b,'ok')






% plot onset response spikes vs pupil
winSz   = 1; 
halfWin = winSz/2;
stepSz  = 0.25;
steps   = -2:stepSz:2.5;

% images
subplot(1,2,1)
a = picPup_ppc;
b = ppcSpkPic;
mn = nan(size(steps)); er = nan(size(steps));
for k = 1:numel(steps)
    win = [steps(k)-halfWin steps(k)+halfWin];
    g = find(a>=win(1)&a<=win(2));
    if numel(g) < 5; continue; end
    mn(k) = nanmean(b(g));
    er(k) = nanstd(b(g))/sqrt(numel(g));
end
errorbar(steps,mn,er,'b')
hold on

a = picPup_pul;
b = pulSpkPic;
mn = nan(size(steps)); er = nan(size(steps));
for k = 1:numel(steps)
    win = [steps(k)-halfWin steps(k)+halfWin];
    g = find(a>=win(1)&a<=win(2));
    if numel(g) < 5; continue; end
    mn(k) = nanmean(b(g));
    er(k) = nanstd(b(g))/sqrt(numel(g));
end
errorbar(steps,mn,er,'r')
hold on
xlim([-2.1 2.7])
% ylim([-0.5 0.5])
xlabel('Pupil diameter at stim onset (z-score)')
ylabel('Z-scored firing rate (first 500ms post-stim)')
title('Images')
legend('PPC','LP/Pulvinar')

% videos
subplot(1,2,2)
a = vidPup_ppc;
b = ppcSpkVid;
mn = nan(size(steps)); er = nan(size(steps));
for k = 1:numel(steps)
    win = [steps(k)-halfWin steps(k)+halfWin];
    g = find(a>=win(1)&a<=win(2));
    if numel(g) < 5; continue; end
    mn(k) = nanmean(b(g));
    er(k) = nanstd(b(g))/sqrt(numel(g));
end
errorbar(steps,mn,er,'b')
hold on

a = vidPup_pul;
b = pulSpkVid;
mn = nan(size(steps)); er = nan(size(steps));
for k = 1:numel(steps)
    win = [steps(k)-halfWin steps(k)+halfWin];
    g = find(a>=win(1)&a<=win(2));
    if numel(g) < 5; continue; end
    mn(k) = nanmean(b(g));
    er(k) = nanstd(b(g))/sqrt(numel(g));
end
errorbar(steps,mn,er,'r')
hold on
xlim([-2.1 2.7])
% ylim([-0.5 0.5])
xlabel('Pupil diameter at stim onset (z-score)')
ylabel('Z-scored firing rate (first 500ms post-stim)')
title('Videos')
legend('PPC','LP/Pulvinar')





% average across units per session
x = reshape(vidPup_ppc,[32 count]);






mn = mean(vidPup);
er = std(vidPup)/size(vidPup,1);

plot(tvec,mn+er,'r--')
plot(tvec,mn-er,'r--')

xlabel('Time (s)')
ylabel('z-scored pupil diameter')
set(gca,'TickDir','out')


figure
subplot(1,2,1)
scatter(picOnPup,picLat)
[r,p] = corrcoef(picOnPup,picLat,'rows','complete');
xlabel('Pupil at image onset (z-score)')
ylabel('First saccade latency')
title(['r = ' num2str(r(1,2)) ' p = ' num2str(p(1,2))])

subplot(1,2,2)
scatter(vidOnPup,vidLat)
[r,p] = corrcoef(vidOnPup,vidLat,'rows','complete');
xlabel('Pupil at video onset (z-score)')
ylabel('First saccade latency')
title(['r = ' num2str(r(1,2)) ' p = ' num2str(p(1,2))])


figure

winSz   = 1; 
halfWin = winSz/2;
stepSz  = 0.25;
steps   = -2.5:stepSz:2.5;

% images
a = picOnPup;
b = picLat;
mn = nan(size(steps)); er = nan(size(steps));
for k = 1:numel(steps)
    win = [steps(k)-halfWin steps(k)+halfWin];
    g = find(a>=win(1)&a<=win(2));
    if numel(g) < 5; continue; end
    mn(k) = nanmean(b(g));
    er(k) = nanstd(b(g))/sqrt(numel(g));
end
errorbar(steps,mn,er,'b')
hold on

a = vidOnPup;
b = vidLat;
mn = nan(size(steps)); er = nan(size(steps));
for k = 1:numel(steps)
    win = [steps(k)-halfWin steps(k)+halfWin];
    g = find(a>=win(1)&a<=win(2));
    if numel(g) < 5; continue; end
    mn(k) = nanmean(b(g));
    er(k) = nanstd(b(g))/sqrt(numel(g));
end
errorbar(steps,mn,er,'r')
hold on

set(gca,'TickDir','out')
legend('Images','Video')
xlabel('Pupil diameter at trial onset (z-score)')
ylabel('First saccade latency (s)')
xlim([-2.7 2.7])    
ylim([0 4.5])




