function  [sacTime,cazi,cele,cpup] = is_detectSaccadesEphysFun_V3(PreprocessDir)
% This function loads in the "ADC" (analogue-to-digital) data from the
% INTAN system. I have the ADC inputs currently set up with the first 3
% channels being ISCAN voltage inputs (eye horizontal position, eye 
% vertical position, pupil diameter). 
% This code automatically detects saccades, and removes blinks from pupil
% diameter traces. 
% I.S. 2017
% AH 2018 V3 switched order of output, added condition to handle 2 eyes,
% i.e. 6 adc channel computation; added recPath input, added figures to
% save into the pupPreprocessDir for inspection
% 10/23/2018 added 30Hz lowpass filter to adc data and then downsample to 1000Hz, 
% change high pass filter that detects blinks from 50Hz to 30Hz 

display(['Detecting saccades ' PreprocessDir])
pupPreprocessDir = [PreprocessDir 'pupilFigures/'];
if ~exist(join([pupPreprocessDir]),'dir'); mkdir(join([pupPreprocessDir]));end % to save figures for each session

lfpFs = is_load([PreprocessDir 'lfp/lfpMat'],'lfpFs');
load([PreprocessDir 'adc_data']);

if size(adc_data,1) < 3 ; display('Incorrect INTAN data'); return; end

%% Filter and resample eye data at the same sample Fs as the lfp
% define a low pass filter at 30Hz, and down sample from 30000 to 1000 Hz
downFac = 30;         % <--------factor to downsample raw signal for LFP
lpFreq  = 30;         % lowpass frequency
ford    = 3;          % filter order
lfpFs   = Fs/downFac; % LFP sample rate
nyqFreq = Fs/2;       % Nyquest frequency of raw signal
[a,b]   = butter(ford,lpFreq/nyqFreq,'low'); % lowpass filter 
[dr,nr] =  rat(Fs/(Fs/downFac)); % Ratio for downsampling
for iChn = 1:size(adc_data,1)
    temp = adc_data(iChn,:);
    filted = filtfilt(a,b,temp);
    dsampled = resample(filted,nr,dr); % downsample data
    adc_data_lfpFs(iChn,:) = dsampled;
end

save([PreprocessDir 'adc_data_lfpFs'],'adc_data_lfpFs','lfpFs');


%% getting eye data
azi   = adc_data_lfpFs(1,:);
ele   = adc_data_lfpFs(2,:);
pup   = adc_data_lfpFs(3,:);

doPlot = 1; 
plotWin = [1,2]; 
plotT = plotWin(1)*60:1/lfpFs:plotWin(2)*60; %in sec
tMask = plotWin(1)*60*lfpFs:plotWin(2)*60*lfpFs;

tvec   = (1:numel(azi))/lfpFs; timeTmp = tvec;
xLim = [tvec(1),tvec(end)];

if doPlot == 1
    fig = figure('name','rEye raw');
    subplot(311);plot(tvec, azi);xlim(xLim); ylabel('azi'); title('raw data'); %1-5min data
    subplot(312);plot(tvec, ele);xlim(xLim); ylabel('ele');
    subplot(313);plot(tvec, pup);xlim(xLim); ylabel('pup');xlabel('Time [sec]');
    saveas(fig, [pupPreprocessDir 'rEye_raw.png']);
end


%% We first need to detect blinks by looking for sharp changes
% hp filter still works better even with data that is already lp filtered at 30Hz
[b,a] = butter(2,30/(lfpFs/2),'high'); %50Hz high pass filter
azZ   = abs(zscore(filtfilt(b,a,azi)));
elZ   = abs(zscore(filtfilt(b,a,ele)));
puZ   = abs(zscore(filtfilt(b,a,pup)));

% % if not hp filter
% azZ   = abs(zscore(azi));
% elZ   = abs(zscore(ele));
% puZ   = abs(zscore(pup));


eyeNames = {'rEye'};
if doPlot == 1
    fig = figure(); 
    yLim = [0,20];
    subplot(311);plot(tvec, azZ);hline(5,'r-');xlim(xLim); ylabel('azi');ylim(yLim);title('z-scored rEye to detect blinks');%1-5min data
    subplot(312);plot(tvec, elZ);hline(5,'r-');xlim(xLim); ylabel('ele');ylim(yLim);
    subplot(313);plot(tvec, puZ);hline(5,'r-');xlim(xLim); ylabel('pup');ylim(yLim);xlabel('Time [sec]');
    saveas(fig, [pupPreprocessDir 'rEye_z-scored.png']);
end


thresh = 0.02*lfpFs; % STD threshold to detect blinks 20ms-40ms
interpSamps             = zeros(size(azZ)); % to merge time points above threshold from az, el and pu
interpSamps(azZ>thresh) = 1;
interpSamps(elZ>thresh) = 1;
interpSamps(puZ>thresh) = 1;
% smooth detected samples with a 2 second window
t      = smooth(interpSamps,round(lfpFs)*2);
aziTmp = azi;
eleTmp = ele;
pupTmp = pup;

delSamps = (t>0); % samples to reject
% delete samples that surpass threshold
aziTmp(delSamps)  = [];
eleTmp(delSamps)  = [];
pupTmp(delSamps)  = [];
timeTmp(delSamps) = [];
% Interpolate out rejected samples
aziClean = interp1(timeTmp,aziTmp,tvec,'linear');
eleClean = interp1(timeTmp,eleTmp,tvec,'linear');
pupClean = interp1(timeTmp,pupTmp,tvec,'linear');

cazi = aziClean;
cele = eleClean;
cpup = pupClean;

if doPlot == 1
    fig = figure(); 
    subplot(311);plot(tvec, aziClean);xlim(xLim); ylabel('azi'); title('rEye blink rejected');%1-5min data
    subplot(312);plot(tvec, eleClean);xlim(xLim); ylabel('ele');
    subplot(313);plot(tvec, pupClean);xlim(xLim); ylabel('pup'); xlabel('Time [sec]');
    saveas(fig, [pupPreprocessDir 'rEye_clean.png']);
end


%% Compute a sliding difference
numSamps = 30; % samples to compute difference (human saccade 80-200ms each = 12~5 Hz)
ele_vel  =  cele(numSamps:end) - cele(1:end-numSamps+1);
azi_vel  =  cazi(numSamps:end) - cazi(1:end-numSamps+1);


stdThresh = 4; % define the number of standard deviations to use as threshold to detect saccades
% compute eye velocity
eye_vel = sqrt(ele_vel.^2+azi_vel.^2);
eye_vel_tmp = eye_vel;
eyeThresh = stdThresh * nanstd(eye_vel(60:end-60)); % velocity
% apply threshold
eye_vel(eye_vel<eyeThresh) = 0;
% find the velocity peaks for saccades
[~,sac] = findpeaks(eye_vel);  % find peak times in ms
sac(sac<10) = [];
sac(sac>numel(eye_vel)-10) = [];

artRejTime = 5; % in samples
preSac     = sac-artRejTime;
postSac    = sac+artRejTime;

output      = cell2mat (arrayfun( @(x,y)(x:y), preSac, postSac, 'UniformOutput', false)); % compute all time points in the sac window (concatenated)
artRejCheck = eye_vel_tmp(output) < eyeThresh; % determine which samples fall under velocity threshold (1's)
artRejInd   = reshape(artRejCheck,[artRejTime*2+1 numel(preSac)]); % reshape into saccades vs sample
artRejSum   = sum(artRejInd,1); % add up values for each saccade
artRejSac   = find(artRejSum > 0); % tag "saccades" that don't have wide enough velocity envelope

sac(artRejSac) = [];

% eliminate saccades that occur within 100ms of each other
isi           = diff(sac);
rejInd        = find(isi<(0.1*lfpFs));
sac(rejInd+1) = [];

sacSamp = {sac};
sacTime = {sac/lfpFs};

if doPlot == 1
    fig = figure('name','saccade_R','Position',[10 50 880 540]); 

    %sacPlot = sac(sac>=plotWin(1)*60*lfpFs & sac<=plotWin(2)*60*lfpFs);
    subplot(411);plot(tvec(1:length(eye_vel_tmp)), eye_vel_tmp);hline(eyeThresh,'r-');
        title(['raw velocity + ' num2str(stdThresh) 'std threshold']);xlim(xLim);
    subplot(412);plot(tvec(1:length(eye_vel_tmp)), eye_vel);title('saccades above threshold');xlim(xLim);
    subplot(413);plot(tvec, cazi);hold on; plot(sacSamp{1}./lfpFs,ones(size(sacSamp{1}))*1.1,'r*');
        title('deblinked,filtered azimuth + final saccades');xlabel('Time [s]');xlim(xLim);
    subplot(414);plot(tvec, cazi);hold on; plot(sacSamp{1}./lfpFs,ones(size(sacSamp{1}))*1.1,'r*');
        title('zoom in -- deblinked,filtered azimuth + final saccades');xlabel('Time [s]');xlim([plotT(1),plotT(end)]);
    %if ~isempty(sacPlot); vline(sacPlot/lfpFs, 'r-');end
    saveas(fig, [pupPreprocessDir 'rEye_saccade.png']);
    savefig(fig, [pupPreprocessDir 'rEye_saccade.fig'],'compact');
end

%% If there is another eye:
% % Filter and resample eye data at the same sample Fs as the lfp
if size(adc_data,1) >= 6
azi   = adc_data_lfpFs(4,:);
ele   = adc_data_lfpFs(5,:);
pup   = adc_data_lfpFs(6,:);
eyeNames = {'lEye','rEye'};
if doPlot == 1
    fig = figure();
    subplot(311);plot(tvec, azi);xlim(xLim);ylabel('azi');title('raw data'); %1-5min data
    subplot(312);plot(tvec, ele);xlim(xLim);ylabel('ele');
    subplot(313);plot(tvec, pup);xlim(xLim);ylabel('pup'); xlabel('Time [sec]');
    saveas(fig, [pupPreprocessDir 'lEye_raw.png']);
end

% We first need to detect blinks by looking for sharp changes
[b,a] = butter(2,30/(lfpFs/2),'high'); %change to 30Hz highpass filter b/c eye data is lp filtered at 30Hz
azZ = abs(zscore(filtfilt(b,a,azi)));
elZ = abs(zscore(filtfilt(b,a,ele)));
puZ = abs(zscore(filtfilt(b,a,pup)));

if doPlot == 1
    fig = figure(); 
    subplot(311);plot(tvec, azZ);hline(5,'r-');xlim(xLim);ylabel('azi'); title('z-scored');%1-5min data
    subplot(312);plot(tvec, elZ);hline(5,'r-');xlim(xLim);ylabel('ele');
    subplot(313);plot(tvec, puZ);hline(5,'r-');xlim(xLim);ylabel('pup'); xlabel('Time [sec]');
    saveas(fig, [pupPreprocessDir 'lEye_z-scored.png']);
end

thresh = 5; % STD threshold to detect blinks
interpSamps             = zeros(size(elZ));
interpSamps(elZ>thresh) = 1;
interpSamps(azZ>thresh) = 1;
interpSamps(puZ>thresh) = 1;
% smooth detected samples with a 2 second window
t      = smooth(interpSamps,round(lfpFs)*2);
aziTmp = azi;
eleTmp = ele;
pupTmp = pup;
tvec   = (1:numel(azi))/lfpFs; timeTmp = tvec;
delSamps = (t>0); % samples to reject
% delete samples that surpass threshold
aziTmp(delSamps)  = [];
eleTmp(delSamps)  = [];
pupTmp(delSamps)  = [];
timeTmp(delSamps) = [];
% Interpolate out rejected samples
aziClean = interp1(timeTmp,aziTmp,tvec,'linear');
eleClean = interp1(timeTmp,eleTmp,tvec,'linear');
pupClean = interp1(timeTmp,pupTmp,tvec,'linear');

cele = [cele; eleClean];
cazi = [cazi; aziClean];
cpup = [cpup; pupClean];

if doPlot == 1
    fig = figure(); 
    subplot(311);plot(tvec, aziClean);xlim(xLim);ylabel('azi'); title('blink rejected');%1-5min data
    subplot(312);plot(tvec, eleClean);xlim(xLim);ylabel('ele');  
    subplot(313);plot(tvec, pupClean);xlim(xLim);ylabel('pup'); xlabel('Time [sec]');
    saveas(fig, [pupPreprocessDir 'lEye_clean.png']);
end


% Compute a sliding difference
numSamps = 30; % samples to compute difference
ele_vel  =  eleClean(numSamps:end) - eleClean(1:end-numSamps+1);
azi_vel  =  aziClean(numSamps:end) - aziClean(1:end-numSamps+1);


stdThresh = 3; % define the number of standard deviations to use as threshold to detect saccades
% compute eye velocity
eye_vel = sqrt(ele_vel.^2+azi_vel.^2);
eye_vel_tmp = eye_vel;
eyeThresh = stdThresh * nanstd(eye_vel(60:end-60)); % elevation
% apply threshold
eye_vel(eye_vel<eyeThresh) = 0;
% find the velocity peaks for saccades
[~,sac] = findpeaks(eye_vel);
sac(sac<10) = [];
sac(sac>numel(eye_vel)-10) = [];

artRejTime = 5; % in samples
preSac     = sac-artRejTime;
postSac    = sac+artRejTime;

output      = cell2mat (arrayfun( @(x,y)(x:y), preSac, postSac, 'UniformOutput', false)); % compute all samples in the sac window (concatenated)
artRejCheck = eye_vel_tmp(output) < eyeThresh; % determine which samples fall under velocity threshold (1's)
artRejInd   = reshape(artRejCheck,[artRejTime*2+1 numel(preSac)]); % reshape into saccades vs sample
artRejSum   = sum(artRejInd,1); % add up values for each saccade
artRejSac   = find(artRejSum > 0); % tag "saccades" that don't have wide enough velocity envelope

sac(artRejSac) = [];

% eliminate saccades that occur within 100ms of each other
isi           = diff(sac);
rejInd        = find(isi<(0.1*lfpFs));
sac(rejInd+1) = [];

sacSamp(2) = {sac};
sacTime(2) = {sac/lfpFs};

% test saccade detecion
if doPlot == 1
    fig = figure('name','saccade_R','Position',[10 50 880 540]); 
    %sacPlot = sac(sac>=plotWin(1)*60*lfpFs & sac<=plotWin(2)*60*lfpFs);
    subplot(411);plot(tvec(1:length(eye_vel_tmp)), eye_vel_tmp);hline(eyeThresh,'r-');
        title(['raw velocity + ' num2str(stdThresh) 'std threshold']);xlim(xLim);
    subplot(412);plot(tvec(1:length(eye_vel_tmp)), eye_vel);title('saccades above threshold');xlim(xLim);
    subplot(413);plot(tvec, cazi);hold on; plot(sacSamp{2}./lfpFs,ones(size(sacSamp{2}))*1.1,'r*');
        title('deblinked,filtered azimuth + final saccades');xlabel('Time [s]');xlim(xLim);legend('rEye','lEye','Saccade');
    subplot(414);plot(tvec, cazi);hold on; plot(sacSamp{2}./lfpFs,ones(size(sacSamp{2}))*1.1,'r*');
        title('zoom in -- deblinked,filtered azimuth + final saccades');xlabel('Time [s]');xlim([plotT(1),plotT(end)]);
    %if ~isempty(sacPlot); vline(sacPlot/lfpFs, 'r-');end
    saveas(fig, [pupPreprocessDir 'lEye_saccade.png']);
    savefig(fig, [pupPreprocessDir 'lEye_saccade.fig'],'compact');
end

cele = flip(cele,1); % flip so that left eye is on top
cazi = flip(cazi,1);
cpup = flip(cpup,1);
sacSamp = flip(sacSamp,1);
sacTime = flip(sacTime,1);
end

%%
save([PreprocessDir 'saccades'],'sacSamp','sacTime')
save([PreprocessDir 'pupil'],'cele','cazi','cpup','eyeNames')

close all
return


% aziClean = filtfilt(b,a,aziClean);
% eleClean = filtfilt(b,a,eleClean);
% 
% % prepare a normalised bartlett window to convolve with data
% smoothTime = round(0.05*lfpFs); % 50ms window
% barWin     = bartlett(smoothTime);
% barWin     = barWin/sum(barWin);
% % smooth data by convolving with bartlett window
% cele   = conv(eleClean,barWin,'same');
% cazi   = conv(aziClean,barWin,'same');

% Check the two-dimensional distribution of eye velocity vectors
% Hout = hist2d(vertcat(azi_vel,ele_vel),100,100,[-0.2 0.2],[-0.2 0.2]);
% imagesc(log10(Hout)); caxis([0 4.5]); colormap(hot)
% xlabel('azimuth velocity')
% ylabel('elevation velocity')

% check eye position
% Hout = hist2d(vertcat(cazi,cele),100,100,[min(cazi) max(cazi)],[min(cele) max(cele)]);
% imagesc(Hout); colormap(hot)
% set(gca,'YDir','normal')
% xlabel('azimuth (a.u)')
% ylabel('elevation (a.u)')
% title('animal pupil center (in pixels)')


% left over code
% [b,a]   = butter(2,40/nyqFreq,'low'); % define a low pass filter at 60Hz
% 
% 
% % interpolate out pupil noise
% df = diff(adc_data(1,:));
% df(df>-0.02 & df<0.02) = 0;
% [negPeaks,negLoc] = findpeaks(df*-1);
% [posPeaks,posLoc] = findpeaks(df);
% delNeg = (negPeaks>0.1); 
% delPos = (posPeaks>0.1);
% negLoc(delNeg) = [];
% posLoc(delPos) = [];
% 
% % plot(adc_data(1,:)); hold on
% % plot(posLoc,ones(size(posLoc)),'r*')
% % plot(negLoc,ones(size(negLoc)),'b*')
% 
% negVec = zeros(1,size(adc_data,2));
% negVec(negLoc) = 1;
% posVec = zeros(1,size(adc_data,2));
% posVec(posLoc) = 1;
% 
% nwin = [zeros(1,200) ones(1,300)];
% pwin = [ones(1,300) zeros(1,200)];
% %bwin(1,:) = bartlett(400); bwin = [zeros(1,200) bwin(201:end)];
% n = conv(negVec,nwin,'same');
% p = conv(posVec,pwin,'same');
% thresh = (n+p >0);
% tvec = (1:size(adc_data,2))/Fs;
% adc1 = adc_data(1,:); adc1(thresh) = [];
% adc2 = adc_data(1,:); adc2(thresh) = [];
% tt   = tvec; tt(thresh) = [];
% intAdc1 = interp1(tt,adc1,tvec,'linear');
%  
% delSamps = [];
% for iev = 1:numel(negLoc)
%     tmpLoc = posLoc;
%     tmp    = posLoc-negLoc(iev);tmpLoc(tmp<0) = []; tmp(tmp<0) = [];
%     ind    = max(tmpLoc(tmp == min(tmp)));
%     if ind < negLoc(iev); continue; end
%     if numel(negLoc(iev):ind) > 30000; continue; end
%     delSamps = [delSamps (negLoc(iev)-50):(ind+50)];
% end
% % plot(adc_data(1,:)); hold on
% % plot(delSamps,ones(size(delSamps)),'r*')
% 
% tvec = (1:size(adc_data,2))/Fs;
% adc1 = adc_data(1,:); adc1(delSamps) = [];
% adc2 = adc_data(2,:); adc2(delSamps) = [];
% adc3 = adc_data(3,:); adc3(delSamps) = [];
% tt   = tvec; tt(delSamps) = [];
% intAdc1 = interp1(tt,adc1,tvec,'linear');
% intAdc2 = interp1(tt,adc2,tvec,'linear');
% intAdc3 = interp1(tt,adc3,tvec,'linear');
% 
% % plot(adc_data(1,:)); hold on
% % plot(intAdc1,'r')
% 
% azi = resample(intAdc1,nr,dr);
% ele = resample(intAdc2,nr,dr);
% pup = resample(intAdc3,nr,dr);




% % compute the inter-saccade interval
% isi        = diff(sac);
% isi(isi<0.05*lfpFs) = []; % reject if isi is less than 6 samples
% histVec = (0.05:0.05:1)*lfpFs;
% bn = histc(isi,histVec)./numel(isi);
% bar(histVec./lfpFs,bn);
% xlabel('inter saccade interval')
% ylabel('probability')
% xlim([0 1])

% % compute a sliding average 
% numSamps = 10;
% convVec  = ones(1,numSamps)/numSamps;
% conv_ele = conv(aziClean,convVec,'same');
% conv_azi = conv(aziClean,convVec,'same');
% % compute velocity on non overlapping data segments
% halfSeg  = round(numSamps/2);
% ele_vel  = conv_ele(halfSeg:end) - cele(1:end-halfSeg+1);
% azi_vel  = conv_azi(halfSeg:end) - cazi(1:end-halfSeg+1);
% % the above convolution and subtraction should be equivalent to computing
% % the sliding average based on the values (numSamps) either side of each
% % data point. 


