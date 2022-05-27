

animals = {'0114','0116','0124','0125'};
pathDir = 'E:\';

bs = 1:4000;

count = 0;
for ianimal = 1:numel(animals)
    animalCode    = animals{ianimal};
    socialDarkDir = dir([pathDir animalCode '\' animalCode '_imagesVideo*']);
    
    for irec = 1:numel(socialDarkDir);
        count = count + 1;
        recName = socialDarkDir(irec).name;
        recPath = [pathDir animalCode '\' recName '\'];
        display(recName)
        
        try
            load([recPath 'analysis\imagesVideoSpectra'])
        catch
            continue
        end
        
        powDiff_vid = nan(size(powVid_ppc,2),size(powVid_ppc,1),size(powVid_ppc,3));
        powDiff_pic = nan(size(powVid_ppc,2),size(powVid_ppc,1),size(powVid_ppc,3));
        for ichan = 1:32
            tmp      = squeeze(powPic_ppc(:,ichan,:));
            baseline = mean(tmp(:,bs),2);
            bsMat    = repmat(baseline,1,size(tmp,2));
            powDiff_pic(ichan,:,:)  = ((tmp-bsMat)./bsMat)*100;
            
            tmp      = squeeze(powVid_ppc(:,ichan,:));
            baseline = mean(tmp(:,bs),2);
            bsMat    = repmat(baseline,1,size(tmp,2));
            powDiff_vid(ichan,:,:)  = ((tmp-bsMat)./bsMat)*100;            
        end
        ppcPow_pic(count,:,:) = squeeze(nanmedian(powDiff_pic));
        ppcPow_vid(count,:,:) = squeeze(nanmedian(powDiff_vid));
        
        powDiff_vid = nan(size(powVid_pul,2),size(powVid_pul,1),size(powVid_pul,3));
        powDiff_pic = nan(size(powVid_pul,2),size(powVid_pul,1),size(powVid_pul,3));
        for ichan = 1:size(powVid_pul,2)
            tmp      = squeeze(powPic_pul(:,ichan,:));
            baseline = mean(tmp(:,bs),2);
            bsMat    = repmat(baseline,1,size(tmp,2));
            powDiff_pic(ichan,:,:)  = ((tmp-bsMat)./bsMat)*100;
            
            tmp      = squeeze(powVid_pul(:,ichan,:));
            baseline = mean(tmp(:,bs),2);
            bsMat    = repmat(baseline,1,size(tmp,2));
            powDiff_vid(ichan,:,:)  = ((tmp-bsMat)./bsMat)*100;            
        end
        pulPow_pic(count,:,:) = squeeze(nanmedian(powDiff_pic));
        pulPow_vid(count,:,:) = squeeze(nanmedian(powDiff_vid));
        
        plv_vid(count,:,:) = avVidPLV;
        plv_pic(count,:,:) = avPicPLV;
    end
end

ppcPow_vid(5:6,:,:) = nan;
pulPow_vid(5:6,:,:) = nan;
ppcPow_pic(5:6,:,:) = nan;
pulPow_pic(5:6,:,:) = nan; pulPow_pic(13,:,:) = nan;

tvec = (-5000:15000)/1000;
lowFreq      = 0.5;
highFreq     = 128;
numFreqs     = 80;
foi          = logspace(log10(lowFreq),log10(highFreq),numFreqs);
fois = [0.5 1 2 4 8 16 32 64 128];
for fi = 1:numel(fois)
    [bi,bb] = sort(abs(foi-fois(fi)));
    tickLoc(fi) = bb(1);
end

mnPPC = squeeze(nanmean(ppcPow_pic));
mnPul = squeeze(nanmean(pulPow_pic));
figure
subplot(1,2,1)
imagesc(tvec,1:numel(foi),mnPPC);
set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',{'0.5','1','2','4','8','16','32','64','128'});
caxis([-60 60])
colormap(jet)
h = colorbar;
ylabel(h,'Relative % power change')
title('PPC Image')
xlabel('Time (s)')
ylabel('Frequency (Hz)')

subplot(1,2,2)
imagesc(tvec,1:numel(foi),mnPul);
set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',{'0.5','1','2','4','8','16','32','64','128'});
caxis([-60 60])
colormap(jet)
h = colorbar;
ylabel(h,'Relative % power change')
title('LP/Pulvinar Image')
xlabel('Time (s)')
ylabel('Frequency (Hz)')


mnPPC = squeeze(nanmean(ppcPow_vid));
mnPul = squeeze(nanmean(pulPow_vid));
figure
subplot(1,2,1)
imagesc(tvec,1:numel(foi),mnPPC);
set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',{'0.5','1','2','4','8','16','32','64','128'});
caxis([-60 60])
colormap(jet)
h = colorbar;
ylabel(h,'Relative % power change')
title('PPC Video')
xlabel('Time (s)')
ylabel('Frequency (Hz)')

subplot(1,2,2)
imagesc(tvec,1:numel(foi),mnPul);
set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',{'0.5','1','2','4','8','16','32','64','128'});
caxis([-60 60])
colormap(jet)
h = colorbar;
ylabel(h,'Relative % power change')
title('LP/Pulvinar Video')
xlabel('Time (s)')
ylabel('Frequency (Hz)')





mnPic = squeeze(nanmean(plv_pic));
mnVid = squeeze(nanmean(plv_vid));

figure
imagesc(tvec,1:numel(foi),mnPic);
set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',{'0.5','1','2','4','8','16','32','64','128'});
caxis([0.15 0.35])
colormap(awesomemap2)
h = colorbar;
ylabel(h,'thalamo-cortical PLV')
title('PPC Video')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
h = colorbar;
ylabel(h,'PLV')

figure
imagesc(tvec,1:numel(foi),mnVid);
set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',{'0.5','1','2','4','8','16','32','64','128'});
caxis([0.15 0.35])
colormap(awesomemap2)
h = colorbar;
ylabel(h,'thalamo-cortical PLV')
title('LP/Pulvinar Video')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
h = colorbar;
ylabel(h,'PLV')

    