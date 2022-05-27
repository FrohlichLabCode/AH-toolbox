
animals = {'0114','0116','0124','0125'};
pathDir      = 'D:\';
lowFreq      = 0.5;
highFreq     = 128;
numFreqs     = 80;
foi          = logspace(log10(lowFreq),log10(highFreq),numFreqs);

nsac = 50;
nrep = 40;
randSamps = 1000;


count = 0;
for ianimal = 1:numel(animals)
    animalCode = animals{ianimal};
    socialDarkDir = dir([pathDir animalCode '\' animalCode '_socialDark*']);
    
    for irec = 1:numel(socialDarkDir)
        recName = socialDarkDir(irec).name;
        recPath = [pathDir animalCode '\' recName '\'];
        display(recName);
        h = keepRecSaccades(recName);
        if h == 0; continue; end
        count = count + 1;
        
        load([recPath 'lfp\lfpMat'])
        load([recPath 'saccades'])
        if numel(sacTime) < nsac; continue; end
        
        keepChans = keepPulChans(animalCode);
        if size(lfpMat,1) >49; pulChans = keepChans + 40; else pulChans = keepChans + 32; end
        ppcChans = 1:32;
        
        ppcSig = median(lfpMat(ppcChans,:));
        pulSig = median(lfpMat(pulChans,:));
        
        wavs = is_makeWavelet(foi,lfpFs);
        win  = round([-10 10]*lfpFs);
        
        %ppcSpec = nan(numel(foi),diff(win)+1);
        %pulSpec = nan(numel(foi),diff(win)+1);
        for f = 1:numel(foi)
            C_ppc = conv(ppcSig,wavs{f},'same');
            C_pul = conv(pulSig,wavs{f},'same');
            phaseDiff = angle(C_ppc) - angle(C_pul);
            
            ppcZ = zscore(abs(C_ppc));
            pulZ = zscore(abs(C_pul));
            
            % loop through each saccade and get spectral data
            sta = nan(length(sacSamp),diff(win)+1);
            sta_ppc = nan(length(sacSamp),diff(win)+1);
            sta_pul = nan(length(sacSamp),diff(win)+1);
            for isac = 1:length(sacSamp)
                if sacSamp(isac) < win(2); continue; end
                if sacSamp(isac) > size(lfpMat,2) - win(2); continue; end
                sta(isac,:) = phaseDiff(:,sacSamp(isac)+win(1):sacSamp(isac)+win(2));
                sta_ppc(isac,:) = ppcZ(:,sacSamp(isac)+win(1):sacSamp(isac)+win(2));
                sta_pul(isac,:) = pulZ(:,sacSamp(isac)+win(1):sacSamp(isac)+win(2));
            end
            
            sacSpec_ppc(f,:) = nanmean(sta_ppc);
            sacSpec_pul(f,:) = nanmean(sta_pul);
            
            ppcSpec(count,f,:) = nanmean(sta_ppc);
            pulSpec(count,f,:) = nanmean(sta_pul);
            
            for irep = 1:nrep
                rp = randperm(size(sta,1));
                tmp_plv(irep,:) = abs(nanmean(exp(1i*sta(rp(1:nsac),:))));
            end
            
            plv(count,f,:) = nanmean(tmp_plv);
            
            % now get random data (for stats later on..)
            recLength   = numel(ppcSig) - 2000;
            plv_rand    = nan(1,randSamps);
            ppcPow_rand = nan(1,randSamps);
            pulPow_rand = nan(1,randSamps);
            for irep = 1:randSamps
                rs = round(rand(1,nsac) * recLength) + 1000;
                plv_rand(irep) = abs(nanmean(exp(1i*phaseDiff(rs))));
                ppcPow_rand(irep) = mean(ppcZ(rs));
                pulPow_rand(irep) = mean(pulZ(rs));
            end
            ppcPowMean(count,f,:) = nanmean(ppcPow_rand);
            ppcPowSTD(count,f,: ) = nanstd(ppcPow_rand);
            pulPowMean(count,f,:) = nanmean(pulPow_rand);
            pulPowSTD(count,f,: ) = nanstd(pulPow_rand);
            plvMean(count,f,:)    = nanmean(plv_rand);
            plvSTD(count,f,:)     = nanstd(plv_rand);
        end
        save([recPath 'analysis\saccadeSpec'],'sacSpec_ppc','sacSpec_pul')
    end
end
        
mn = squeeze(nanmean(plv));
tvec = (win(1):win(2))/lfpFs;
fois = [0.5 1 2 4 8 16 32 64 128];
for fi = 1:numel(fois)
    [bi,bb] = sort(abs(foi-fois(fi)));
    tickLoc(fi) = bb(1);
end

figure
imagesc(tvec,1:numel(foi),mn)
set(gca,'TickDir','out','YDir','normal','YTick',tickLoc,'YTickLabel',{'0.5','1','2','4','8','16','32','64','128'});
colormap(awesomeMap)
caxis([0.1 0.45])

figure
imagesc(tvec,1:numel(foi),mn)
set(gca,'TickDir','out','YDir','normal','YTick',tickLoc,'YTickLabel',{'0.5','1','2','4','8','16','32','64','128'});
colormap(awesomeMap)
caxis([0.10 0.23])

plvNorm = nan(size(plv));
ppcNorm = nan(size(ppcSpec));
pulNorm = nan(size(pulSpec));
for irec = 1:size(plv,1)
    mnMat = repmat(plvMean(irec,:)',1,size(plv,3));
    stMat = repmat(plvSTD(irec,:)',1,size(plv,3));
    plvNorm(irec,:,:) = (squeeze(plv(irec,:,:)) - mnMat)./stMat;

    mnMat = repmat(ppcPowMean(irec,:)',1,size(ppcPowMean,3));
    stMat = repmat(ppcPowSTD(irec,:)',1,size(ppcPowMean,3));
    ppcNorm(irec,:,:) = (squeeze(ppcSpec(irec,:,:)) - mnMat)./stMat;
    
    mnMat = repmat(pulPowMean(irec,:)',1,size(pulPowMean,3));
    stMat = repmat(pulPowSTD(irec,:)',1,size(pulPowMean,3));
    pulNorm(irec,:,:) = (squeeze(pulSpec(irec,:,:)) - mnMat)./stMat;
end

figure
subplot(1,2,1)
imagesc(tvec,1:numel(foi),squeeze(sum(plvNorm > 2))/size(plv,1) * 100)
set(gca,'TickDir','out','YDir','normal','YTick',tickLoc,'YTickLabel',{'0.5','1','2','4','8','16','32','64','128'});
colormap(awesomeMap)
caxis([0 50])
xlabel('time to saccade')
ylabel('Frequency')
title('PLV significant increase')

subplot(1,2,2)
imagesc(tvec,1:numel(foi),squeeze(sum(plvNorm < -2))/size(plv,1) * 100)
set(gca,'TickDir','out','YDir','normal','YTick',tickLoc,'YTickLabel',{'0.5','1','2','4','8','16','32','64','128'});
colormap(awesomeMap)
caxis([0 50])
xlabel('time to saccade')
ylabel('Frequency')
title('PV significant decrease')

figure
subplot(2,2,1)
imagesc(tvec,1:numel(foi),squeeze(sum(ppcNorm > 2))/size(plv,1) * 100)
set(gca,'TickDir','out','YDir','normal','YTick',tickLoc,'YTickLabel',{'0.5','1','2','4','8','16','32','64','128'});
colormap(awesomeMap)
caxis([0 80])
title('PPC sig increase')
xlabel('time to saccade')
ylabel('Frequency')

subplot(2,2,2)
imagesc(tvec,1:numel(foi),squeeze(sum(ppcNorm < -2))/size(plv,1) * 100)
set(gca,'TickDir','out','YDir','normal','YTick',tickLoc,'YTickLabel',{'0.5','1','2','4','8','16','32','64','128'});
colormap(awesomeMap)
caxis([0 80])
title('PPC sig decrease')
xlabel('time to saccade')
ylabel('Frequency')

subplot(2,2,3)
imagesc(tvec,1:numel(foi),squeeze(sum(pulNorm > 2))/size(plv,1) * 100)
set(gca,'TickDir','out','YDir','normal','YTick',tickLoc,'YTickLabel',{'0.5','1','2','4','8','16','32','64','128'});
colormap(awesomeMap)
caxis([0 80])
title('LP/Pulvinar sig increase')
xlabel('time to saccade')
ylabel('Frequency')

subplot(2,2,4)
imagesc(tvec,1:numel(foi),squeeze(sum(pulNorm < -2))/size(plv,1) * 100)
set(gca,'TickDir','out','YDir','normal','YTick',tickLoc,'YTickLabel',{'0.5','1','2','4','8','16','32','64','128'});
colormap(awesomeMap)
caxis([0 80])
title('LP/Pulvinar sig decrease')
xlabel('time to saccade')
ylabel('Frequency')


bl = (foi<8);
ab = (foi>8);



subplot(2,1,2)
imagesc(tvec,20:39,mn(bl,:))
set(gca,'TickDir','out','YDir','normal','YTick',tickLoc(3:4),'YTickLabel',{'2','4'});
colormap(awesomeMap)
caxis([0.15 0.40])
ylim([20 39])
colorbar

subplot(2,1,1)
imagesc(tvec,40:64,mn(ab,:))
set(gca,'TickDir','out','YDir','normal','YTick',tickLoc(5:7),'YTickLabel',{'8','16','32'});
colormap(awesomeMap)
caxis([0.12 0.25])
ylim([40 64])
colorbar
