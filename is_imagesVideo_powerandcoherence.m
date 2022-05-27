
animals = {'0114','0116','0124','0125'};
pathDir       = 'D:\';
lowFreq      = 0.5;
highFreq     = 128;
numFreqs     = 80;
foi          = logspace(log10(lowFreq),log10(highFreq),numFreqs);
logF         = logspace(log10(lowFreq),log10(highFreq),1000);
newFs        = 200;

count = 0;

for ianimal = 1:numel(animals)
    animalCode = animals{ianimal};
    socialDarkDir = dir([pathDir animalCode '\' animalCode '_imagesVideo*']);
    
    for irec = 1:numel(socialDarkDir);
        
        recName = socialDarkDir(irec).name;
        recPath = [pathDir animalCode '\' recName '\'];
        load([recPath 'triggerData'])
        load([recPath 'saccades'])
        
        try
            % read in the log file
            fileName = ['D:\Presentation_LogFiles\' recName(1:end-7) '.log'];
            fileID   = fopen(fileName);
            formatSpec = '%f %s %s %f %f %f %f %f %f %f %s %f';
            LogFile = textscan(fileID,formatSpec,'HeaderLines',5,'Delimiter', '\t');
            fclose(fileID);
        catch; continue; end
        
        display(recName)
        count = count + 1;
        load([recPath 'lfp\lfpMat'])
        keepChans = keepPulChans(animalCode);
        if size(lfpMat,1) >48; pulChans = keepChans + 40; else pulChans = keepChans + 32; end
        ppcChans = 1:32;
        
        % detect images
        cellPic = cellfun(@regexp,LogFile(2),{'Picture'},'UniformOutput',false);
        picInd  = find(cellfun(@isempty,cellPic{:}) == 0);
        % detect videos
        cellVid = cellfun(@regexp,LogFile(2),{'Video'},'UniformOutput',false);
        vidInd  = find(cellfun(@isempty,cellVid{:}) == 0);
        
        % detect stim onsets
        onset  = find(diff(triggerData(1,:))==1);
        offset = find(diff(triggerData(2,:))==1);
        
        picOnset  = (onset(picInd)/Fs);
        vidOnset  = (onset(vidInd)/Fs);
        %         picOffset = (offset(picInd)/Fs);
        %         vidOffset = (offset(vidInd)/Fs);
        
        window   = 1024;
        noverlap = 512;
        
        clear ppcMat
        for ichan = 1:numel(ppcChans)
            [s,f,t] = spectrogram(lfpMat(ppcChans(ichan),:),window,noverlap,foi,lfpFs);
            ppcMat(ichan,:,:) = s;
        end
        
        clear pulMat
        for ichan = 1:numel(pulChans)
            [s,f,t] = spectrogram(lfpMat(pulChans(ichan),:),window,noverlap,foi,lfpFs);
            pulMat(ichan,:,:) = s;
        end
        
        fs = lfpFs/noverlap;
        gs = 1:10*fs;
        
        % compute pic samps 
        picSamps = [];
        for iev = 1:numel(picOnset)
            picSamps = [picSamps round(fs*picOnset(iev))+gs];
        end
        picSamps(picSamps>numel(t)) = [];
        
        % compute vid samps 
        vidSamps = [];
        for iev = 1:numel(vidOnset)
            vidSamps = [vidSamps round(fs*vidOnset(iev))+gs];
        end        
        vidSamps(vidSamps>numel(t)) = [];
        
        restSamps = 1:numel(t);
        restSamps(ismember(restSamps,[picSamps vidSamps])) = [];
        
        ppcVidPow(count,:) = squeeze(mean(nanmean(abs(ppcMat(:,:,vidSamps)).^2,3)));
        ppcPicPow(count,:) = squeeze(mean(nanmean(abs(ppcMat(:,:,picSamps)).^2,3)));
        ppcGraPow(count,:) = squeeze(mean(nanmean(abs(ppcMat(:,:,restSamps)).^2,3)));
        pulVidPow(count,:) = squeeze(mean(nanmean(abs(pulMat(:,:,vidSamps)).^2,3)));
        pulPicPow(count,:) = squeeze(mean(nanmean(abs(pulMat(:,:,picSamps)).^2,3)));
        pulGraPow(count,:) = squeeze(mean(nanmean(abs(pulMat(:,:,restSamps)).^2,3)));
        
        % coherence video
        for ichan = 1:numel(ppcChans)
            xdat = squeeze(ppcMat(ichan,:,:));
            for jchan = 1:numel(pulChans)
                ydat = squeeze(pulMat(jchan,:,:));
                
                Sxy(:,1) = nanmean(xdat(:,vidSamps).*conj(ydat(:,vidSamps)),2);
                Sxx(:,1) = nanmean(xdat(:,vidSamps).*conj(xdat(:,vidSamps)),2);
                Syy(:,1) = nanmean(ydat(:,vidSamps).*conj(ydat(:,vidSamps)),2);
                Cy(ichan,jchan,:) = Sxy./sqrt(Sxx.*Syy);
            end
        end
        for f = 1:numFreqs
            tmp = squeeze(abs(Cy(:,:,f)));
            vidCy(count,f) = nanmean(tmp(:));
        end
     
        % coherence images
        for ichan = 1:numel(ppcChans)
            xdat = squeeze(ppcMat(ichan,:,:));
            for jchan = 1:numel(pulChans)
                ydat = squeeze(pulMat(jchan,:,:));
                
                Sxy(:,1) = nanmean(xdat(:,picSamps).*conj(ydat(:,picSamps)),2);
                Sxx(:,1) = nanmean(xdat(:,picSamps).*conj(xdat(:,picSamps)),2);
                Syy(:,1) = nanmean(ydat(:,picSamps).*conj(ydat(:,picSamps)),2);
                Cy(ichan,jchan,:) = Sxy./sqrt(Sxx.*Syy);
            end
        end
        for f = 1:numFreqs
            tmp = squeeze(abs(Cy(:,:,f)));
            picCy(count,f) = nanmean(tmp(:));
        end        
        
        % coherence gray
        for ichan = 1:numel(ppcChans)
            xdat = squeeze(ppcMat(ichan,:,:));
            for jchan = 1:numel(pulChans)
                ydat = squeeze(pulMat(jchan,:,:));
                
                Sxy(:,1) = nanmean(xdat(:,restSamps).*conj(ydat(:,restSamps)),2);
                Sxx(:,1) = nanmean(xdat(:,restSamps).*conj(xdat(:,restSamps)),2);
                Syy(:,1) = nanmean(ydat(:,restSamps).*conj(ydat(:,restSamps)),2);
                Cy(ichan,jchan,:) = Sxy./sqrt(Sxx.*Syy);
            end
        end
        for f = 1:numFreqs
            tmp = squeeze(abs(Cy(:,:,f)));
            grayCy(count,f) = nanmean(tmp(:));
        end                
    end
end

ticks = [0.5 1 2 4 8 16 32 64 128];

subplot(1,3,1)
loglog(foi,mean(ppcGraPow),'k'); hold on
loglog(foi,mean(ppcPicPow),'m')
loglog(foi,mean(ppcVidPow),'r')
xlabel('Frequency')
ylabel('Power')
title('PPC')
legend('Gray','Image','Video')
xlim([foi(1) foi(end)])
set(gca,'TickDir','out','XTick',ticks)

subplot(1,3,2)
loglog(foi,mean(pulGraPow),'k'); hold on
loglog(foi,mean(pulPicPow),'m')
loglog(foi,mean(pulVidPow),'r')
xlabel('Frequency')
ylabel('Power')
title('LP/Pulvinar')
legend('Gray','Image','Video')
xlim([foi(1) foi(end)])
set(gca,'TickDir','out','XTick',ticks)

subplot(1,3,3)
semilogx(foi,mean(grayCy),'k'); hold on
semilogx(foi,mean(picCy),'m')
semilogx(foi,mean(vidCy),'r')
xlabel('Frequency')
ylabel('Coherence')
legend('Gray','Image','Video')
xlim([foi(1) foi(end)])
set(gca,'TickDir','out','XTick',ticks)





        