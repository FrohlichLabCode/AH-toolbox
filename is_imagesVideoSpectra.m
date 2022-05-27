function is_imagesVideoSpectra(recPath,recName,animalCode)

% define frequency parameters
lowFreq      = 0.5;
highFreq     = 128;
numFreqs     = 80;
foi          = logspace(log10(lowFreq),log10(highFreq),numFreqs);

% load LFP
load([recPath 'lfp\lfpMat'])
load([recPath 'triggerData'])
keepChans = keepPulChans(animalCode);
pulChans = 
if size(lfpMat,1) >48; pulChans = keepChans + 40; else pulChans = keepChans + 32; end
ppcChans = 1:32;

% get wavelets
wavs = is_makeWavelet(foi,lfpFs);

try
    % read in the log file
    fileName = ['E:\Presentation_LogFiles\' recName(1:end-7) '.log'];
    fileID   = fopen(fileName);
    formatSpec = '%f %s %s %f %f %f %f %f %f %f %s %f';
    LogFile = textscan(fileID,formatSpec,'HeaderLines',5,'Delimiter', '\t');
    fclose(fileID);
catch; return; end

% detect images
cellPic = cellfun(@regexp,LogFile(2),{'Picture'},'UniformOutput',false);
picInd  = cellfun(@isempty,cellPic{:}) == 0;
% detect videos
cellVid = cellfun(@regexp,LogFile(2),{'Video'},'UniformOutput',false);
vidInd  = cellfun(@isempty,cellVid{:}) == 0;

% detect stim onsets
onset  = find(diff(triggerData(1,:))==1);
offset = find(diff(triggerData(2,:))==1);

% define window
win = round([-5 15]*lfpFs);

picOnset = round((onset(picInd)/Fs)*lfpFs); picOnset(picOnset < abs(win(1)) | picOnset > size(lfpMat,2)+win(1)) = [];
vidOnset = round((onset(vidInd)/Fs)*lfpFs); vidOnset(vidOnset < abs(win(1)) | vidOnset > size(lfpMat,2)+win(1)) = [];

% initialise matrices
powPic_ppc = nan(length(foi),numel(ppcChans),diff(win)+1); 
powPic_pul = nan(length(foi),numel(pulChans),diff(win)+1); 
powVid_ppc = nan(length(foi),numel(ppcChans),diff(win)+1); 
powVid_pul = nan(length(foi),numel(pulChans),diff(win)+1); 
for f = 1:numel(foi)
    display(['img-video ' recName(1:end-7) ' ' num2str(f) '/' num2str(numFreqs)])
    % PPC first
    C_ppc = conv2(lfpMat(ppcChans,:),wavs{f},'same'); % convole PPC LFP with wavelet
    zpow  = zscore(abs(C_ppc).^2,[],2); % z-score power
    angPPC = angle(C_ppc); clear C_ppc
    
    % Images PPC
    smat = nan(length(picOnset),size(zpow,1),diff(win)+1);
    for iev = 1:numel(picOnset); smat(iev,:,:) = zpow(:,picOnset(iev)+win(1):picOnset(iev)+win(2)); end
    powPic_ppc(f,:,:) = squeeze(nanmean(smat));
    
    % Video PPC
    smat = nan(length(vidOnset),size(zpow,1),diff(win)+1);
    for iev = 1:numel(vidOnset); smat(iev,:,:) = zpow(:,vidOnset(iev)+win(1):vidOnset(iev)+win(2)); end
    powVid_ppc(f,:,:) = squeeze(nanmean(smat));

    % Pulvinar second
    C_pul = conv2(lfpMat(pulChans,:),wavs{f},'same'); % convole PPC LFP with wavelet
    zpow  = zscore(abs(C_pul).^2,[],2); % z-score power
    angPul = angle(C_pul); clear C_pul
    
    % Images Pul
    smat = nan(length(picOnset),size(zpow,1),diff(win)+1);
    for iev = 1:numel(picOnset); smat(iev,:,:) = zpow(:,picOnset(iev)+win(1):picOnset(iev)+win(2)); end
    powPic_pul(f,:,:) = squeeze(nanmean(smat));
    
    % Video Pul
    smat = nan(length(vidOnset),size(zpow,1),diff(win)+1);
    for iev = 1:numel(vidOnset); smat(iev,:,:) = zpow(:,vidOnset(iev)+win(1):vidOnset(iev)+win(2)); end
    powVid_pul(f,:,:) = squeeze(nanmean(smat));
    
    clear zpow
    
    % compute PLV
    plvWin = round([-0.5 0.5]*lfpFs);
    steps  = -5:0.01:15;
    % define image samples to use
    picSamps = cell(1,numel(steps));
    for k = 1:numel(steps);
        ts = [];
        for iev = 1:numel(picOnset)
            ts = [ts picOnset(iev)+plvWin(1)+steps(k)*lfpFs:picOnset(iev)+plvWin(2)+steps(k)*lfpFs];
        end
        picSamps{k} = ts;
    end
    % define video samples to use
    vidSamps = cell(1,numel(steps));
    for k = 1:numel(steps);
        ts = [];
        for iev = 1:numel(vidOnset)
            ts = [ts vidOnset(iev)+plvWin(1)+steps(k)*lfpFs:vidOnset(iev)+plvWin(2)+steps(k)*lfpFs];
        end
        vidSamps{k} = ts;
    end
    
    pmat   = nan(numel(ppcChans)*numel(pulChans),numel(steps));
    vmat   = nan(numel(ppcChans)*numel(pulChans),numel(steps));
    count = 0;
    for ichan = 1:numel(ppcChans)
        display(num2str(ichan))
        for jchan = 1:numel(pulChans)
            count = count + 1;
            for istep = 1:numel(steps)
                phaseDiff = angPPC(ichan,picSamps{istep}) - angPul(jchan,picSamps{istep});
                pmat(count,istep) = abs(nanmean(exp(1i*phaseDiff)));
                phaseDiff = angPPC(ichan,vidSamps{istep}) - angPul(jchan,vidSamps{istep});
                vmat(count,istep) = abs(nanmean(exp(1i*phaseDiff)));                
            end
        end
    end
    avPicPLV(f,:) = nanmean(pmat);
    avVidPLV(f,:) = nanmean(vmat);
    
    save([recPath 'analysis\imagesVideoSpectra'],'avPicPLV','avVidPLV','powPic_ppc','powPic_pul','powVid_ppc','powVid_pul','-v7.3');
end




