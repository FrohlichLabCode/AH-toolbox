
animals = {'0114','0116','0124','0125'};
pathDir       = 'D:\';
lowFreq      = 0.5;
highFreq     = 128;
numFreqs     = 80;
foi          = logspace(log10(lowFreq),log10(highFreq),numFreqs);
logF         = logspace(log10(lowFreq),log10(highFreq),1000);
newFs        = 200;

count = 0;
tPLV_pic = [];
aPLV_pic = [];
sac_pic  = [];
lum_pic  = [];
tPLV_vid = [];
aPLV_vid = [];
sac_vid  = [];
lum_vid  = [];

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
        
        intensityVar = imagesVideo_stimuliLuminance;

        % get the luminance values for each stimulus per condition
        % images
        for iev = 1:numel(picInd)
            picName = LogFile{3}(picInd(iev));
            testPic = cellfun(@regexp,{intensityVar(:,1)},picName,'UniformOutput',false);
            testInd = find(cellfun(@isempty,testPic{:}) == 0);
            picLum(iev) = intensityVar{testInd,2};
        end
        % videos
        for iev = 1:numel(vidInd)
            vidName = LogFile{3}(vidInd(iev));
            testVid = cellfun(@regexp,{intensityVar(:,1)},vidName,'UniformOutput',false);
            testInd = find(cellfun(@isempty,testVid{:}) == 0);
            vidLum(iev) = intensityVar{testInd,2};
        end            
        
        [b,a] = butter(2,100/(lfpFs/2),'low');
        fdat  = filtfilt(b,a,lfpMat')';
        
        idat = resample(median(fdat(ppcChans,:)),newFs,lfpFs);
        jdat = resample(median(fdat(pulChans,:)),newFs,lfpFs);
        
        wavs = is_makeWavelet([4 14],newFs);
        
        itheta = conv(idat,wavs{1},'same');
        jtheta = conv(jdat,wavs{1},'same');
        ialpha = conv(idat,wavs{2},'same');
        jalpha = conv(jdat,wavs{2},'same');
        
        numSac = []; thetaPLV = []; alphaPLV = []; lumLev = [];
        for iev = 1:numel(vidOnset)
            numSac(iev)    = sum(sacTime>vidOnset(iev)&sacTime<vidOnset(iev)+10);
            grabSamps = round((vidOnset(iev))*newFs):round((vidOnset(iev) + 10)*newFs);
            phaseDiff = angle(itheta(grabSamps)) - angle(jtheta(grabSamps));
            thetaPLV(iev) = abs(mean(exp(1i*phaseDiff)));
            phaseDiff = angle(ialpha(grabSamps)) - angle(jalpha(grabSamps));
            alphaPLV(iev) = abs(mean(exp(1i*phaseDiff)));
            lumLev(iev) = vidLum(iev);
        end
        
        tPLV_vid = [tPLV_vid thetaPLV];
        aPLV_vid = [aPLV_vid alphaPLV];
        sac_vid  = [sac_vid numSac];   
        lum_vid  = [lum_vid lumLev]; 
        
        numSac = []; thetaPLV = []; alphaPLV = []; lumLev = [];
        for iev = 1:numel(picOnset)
            numSac(iev)    = sum(sacTime>picOnset(iev)&sacTime<picOnset(iev)+10);
            grabSamps = round((picOnset(iev))*newFs):round((picOnset(iev) + 10)*newFs);
            phaseDiff = angle(itheta(grabSamps)) - angle(jtheta(grabSamps));
            thetaPLV(iev) = abs(mean(exp(1i*phaseDiff)));
            phaseDiff = angle(ialpha(grabSamps)) - angle(jalpha(grabSamps));
            alphaPLV(iev) = abs(mean(exp(1i*phaseDiff)));
            lumLev(iev) = picLum(iev);
        end
        
        tPLV_pic = [tPLV_pic (thetaPLV)];
        aPLV_pic = [aPLV_pic (alphaPLV)];
        sac_pic  = [sac_pic numSac];  
        lum_pic  = [lum_pic lumLev]; 
    end
end
xvec = 0:0.01:1;


figure
subplot(1,2,1)
scatter(tPLV_pic,sac_pic,'ok'); hold on
[r,m,b] = regression(tPLV_pic,sac_pic);
y = xvec*m + b;
plot(xvec,y,'r--')
xlabel('Theta PLV')
ylabel('Number of saccades per trial')
[r,p] = corrcoef(tPLV_pic,sac_pic);
title(['Images  - corr coef = ' num2str(r(1,2)) ' p = ' num2str(p(1,2))])
xlim([0 0.9])
ylim([0 25])
set(gca,'TickDir','out')

subplot(1,2,2)
scatter(aPLV_pic,sac_pic,'ok'); hold on
[r,m,b] = regression(aPLV_pic,sac_pic);
y = xvec*m + b;
plot(xvec,y,'r--')
xlabel('Alpha PLV')
ylabel('Number of saccades per trial')
[r,p] = corrcoef(aPLV_pic,sac_pic);
title(['Images - corr coef = ' num2str(r(1,2)) ' p = ' num2str(p(1,2))])
xlim([0 0.9])
ylim([0 25])
set(gca,'TickDir','out')

figure
subplot(1,2,1)
scatter(tPLV_vid,sac_vid,'ok'); hold on
[r,m,b] = regression(tPLV_vid,sac_vid);
y = xvec*m + b;
plot(xvec,y,'r--')
xlabel('Theta PLV')
ylabel('Number of saccades per trial')
[r,p] = corrcoef(tPLV_vid,sac_vid);
title(['Video - corr coef = ' num2str(r(1,2)) ' p = ' num2str(p(1,2))])
xlim([0 0.9])
ylim([0 25])
set(gca,'TickDir','out')

subplot(1,2,2)
scatter(aPLV_vid,sac_vid,'ok')
[r,m,b] = regression(aPLV_vid,sac_vid); hold on
y = xvec*m + b;
plot(xvec,y,'r--')
xlabel('Alpha PLV')
ylabel('Number of saccades per trial')
[r,p] = corrcoef(aPLV_vid,sac_vid);
title(['Video - corr coef = ' num2str(r(1,2)) ' p = ' num2str(p(1,2))])
xlim([0 0.9])
ylim([0 25])
set(gca,'TickDir','out')

[r,m,b] = regression(sac_pic,tPLV_pic);
plotregression(tPLV_pic,sac_pic)

% combined
a = [tPLV_pic tPLV_vid];
b = [sac_pic sac_vid];

subplot(1,2,1)
scatter(a,b)
xlabel('Theta PLV')
ylabel('Number of saccades per trial')
[r,p] = corrcoef(a,b);
title(['corr coef = ' num2str(r(1,2)) ' p = ' num2str(p(1,2))])

a = [aPLV_pic aPLV_vid];
b = [sac_pic sac_vid];

subplot(1,2,2)
scatter(a,b)
xlabel('Alpha PLV')
ylabel('Number of saccades per trial')
[r,p] = corrcoef(a,b);
title(['corr coef = ' num2str(r(1,2)) ' p = ' num2str(p(1,2))])


winSz   = 0.2; 
halfWin = winSz/2;
stepSz  = 0.05;
steps   = winSz/2:stepSz:1-winSz/2;

% images
subplot(1,2,1)
a = tPLV_pic;
b = sac_pic;
mn = nan(size(steps)); er = nan(size(steps));
for k = 1:numel(steps)
    win = [steps(k)-halfWin steps(k)+halfWin];
    g = find(a>=win(1)&a<=win(2));
    if numel(g) < 20; continue; end
    mn(k) = mean(b(g));
    er(k) = std(b(g))/sqrt(numel(g));
end
errorbar(steps,mn,er,'ob')
hold on
[r,m,b] = regression(a,b);
y = xvec*m + b;
plot(xvec,y,'b--')


a = aPLV_pic;
b = sac_pic;
mn = nan(size(steps)); er = nan(size(steps));
for k = 1:numel(steps)
    win = [steps(k)-halfWin steps(k)+halfWin];
    g = find(a>=win(1)&a<=win(2));
    if numel(g) < 20; continue; end
    mn(k) = mean(b(g));
    er(k) = std(b(g))/sqrt(numel(g));
end
errorbar(steps,mn,er,'or')
hold on
[r,m,b] = regression(a,b);
y = xvec*m + b;
plot(xvec,y,'r--')

set(gca,'TickDir','out')
legend('Theta','Alpha')
xlabel('PLV')
ylabel('Number of saccades per trial')
xlim([0 0.8])    
ylim([0 8])
title('Images')
    
% videos
subplot(1,2,2)
a = tPLV_vid;
b = sac_vid;
mn = nan(size(steps)); er = nan(size(steps));
for k = 1:numel(steps)
    win = [steps(k)-halfWin steps(k)+halfWin];
    g = find(a>=win(1)&a<=win(2));
    if numel(g) < 20; continue; end
    mn(k) = mean(b(g));
    er(k) = std(b(g))/sqrt(numel(g));
end
errorbar(steps,mn,er,'ob')
hold on
[r,m,b] = regression(a,b);
y = xvec*m + b;
plot(xvec,y,'b--')

a = aPLV_vid;
b = sac_vid;
mn = nan(size(steps)); er = nan(size(steps));
for k = 1:numel(steps)
    win = [steps(k)-halfWin steps(k)+halfWin];
    g = find(a>=win(1)&a<=win(2));
    if numel(g) < 20; continue; end
    mn(k) = mean(b(g));
    er(k) = std(b(g))/sqrt(numel(g));
end
errorbar(steps,mn,er,'or')
hold on
[r,m,b] = regression(a,b);
y = xvec*m + b;
plot(xvec,y,'r--')

set(gca,'TickDir','out')
legend('Theta','Alpha')
xlabel('PLV')
ylabel('Number of saccades per trial')
xlim([0 0.8]) 
ylim([0 8])
title('Video')


% sort based on saccades per trial

steps   = 1:15;

% images
subplot(1,2,1)
a = tPLV_pic;
b = sac_pic;
mn = nan(size(steps)); er = nan(size(steps));
for k = 1:numel(steps)
    g = find(b==k);
    if numel(g) < 5; continue; end
    mn(k) = mean(a(g));
    er(k) = std(a(g))/sqrt(numel(g));
end
errorbar(steps,mn,er,'b')
hold on

a = aPLV_pic;
b = sac_pic;
mn = nan(size(steps)); er = nan(size(steps));
for k = 1:numel(steps)
    g = find(b==k);
    if numel(g) < 5; continue; end
    mn(k) = mean(a(g));
    er(k) = std(a(g))/sqrt(numel(g));
end
errorbar(steps,mn,er,'r')
hold on

set(gca,'TickDir','out')
legend('Theta','Alpha')
xlabel('Number of saccades per trial')
ylabel('PLV')
xlim([0 11])    
ylim([0 0.6])
title('Images')
    
% videos
subplot(1,2,2)
a = tPLV_vid;
b = sac_vid;
mn = nan(size(steps)); er = nan(size(steps));
for k = 1:numel(steps)
    g = find(b==k);
    if numel(g) < 5; continue; end
    mn(k) = mean(a(g));
    er(k) = std(a(g))/sqrt(numel(g));
end
errorbar(steps,mn,er,'b')
hold on

a = aPLV_vid;
b = sac_vid;
mn = nan(size(steps)); er = nan(size(steps));
for k = 1:numel(steps)
    g = find(b==k);
    if numel(g) < 5; continue; end
    mn(k) = mean(a(g));
    er(k) = std(a(g))/sqrt(numel(g));
end
errorbar(steps,mn,er,'r')
hold on

set(gca,'TickDir','out')
legend('Theta','Alpha')
xlabel('Number of saccades per trial')
ylabel('PLV')
xlim([0 15]) 
ylim([0 0.6])
title('Video')


% look at relationship between luminance and PLV
figure
subplot(2,2,1)
scatter(lum_pic,tPLV_pic)
[r,p] = corrcoef(lum_pic,tPLV_pic);
xlabel('Image luminance')
ylabel('thalamo-cortical theta PLV')
title(['r = ' num2str(r(1,2)) ' p = ' num2str(p(1,2))])
xlim([34 46])

subplot(2,2,2)
scatter(lum_pic,aPLV_pic)
[r,p] = corrcoef(lum_pic,aPLV_pic);
xlabel('Image luminance')
ylabel('thalamo-cortical alpha PLV')
title(['r = ' num2str(r(1,2)) ' p = ' num2str(p(1,2))])
xlim([34 46])

subplot(2,2,3)
scatter(lum_vid,tPLV_vid)
[r,p] = corrcoef(lum_vid,tPLV_vid);
xlabel('Video luminance')
ylabel('thalamo-cortical theta PLV')
title(['r = ' num2str(r(1,2)) ' p = ' num2str(p(1,2))])
xlim([34 46])

subplot(2,2,4)
scatter(lum_vid,aPLV_vid)
[r,p] = corrcoef(lum_vid,aPLV_vid);
xlabel('Image luminance')
ylabel('thalamo-cortical alpha PLV')
title(['r = ' num2str(r(1,2)) ' p = ' num2str(p(1,2))])
xlim([34 46])

% number of saccades and luminance
figure
subplot(1,2,1)
scatter(lum_pic,sac_pic)
[r,p] = corrcoef(lum_pic,sac_pic);
xlabel('Image luminance')
ylabel('thalamo-cortical theta PLV')
title(['r = ' num2str(r(1,2)) ' p = ' num2str(p(1,2))])
xlim([34 46])

subplot(1,2,2)
scatter(lum_vid,sac_vid)
[r,p] = corrcoef(lum_vid,sac_vid);
xlabel('Image luminance')
ylabel('thalamo-cortical alpha PLV')
title(['r = ' num2str(r(1,2)) ' p = ' num2str(p(1,2))])
xlim([34 46])



figure
subplot(1,2,1)
hist(sac_pic,0:1:40)
xlabel('# saccades')
ylabel('count')
title('Images')
xlim([-0.5 20])
subplot(1,2,2)
hist(sac_vid,0:1:40)
xlabel('# saccades')
ylabel('count')
title('Video')
xlim([-0.5 20])




