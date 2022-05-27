function is_imagesVideoSpikePLV(recPath,recName,animalCode)

% define frequency parameters
lowFreq      = 0.5;
highFreq     = 128;
numFreqs     = 80;
foi          = logspace(log10(lowFreq),log10(highFreq),numFreqs);

% load LFP
is_cleanSpikes(recPath)
load([recPath 'lfp/lfpMat'])
load([recPath 'triggerData'])
keepChans = keepPulChans(animalCode);
if size(lfpMat,1) >48; pulChans = keepChans + 40; else pulChans = keepChans + 32; end
ppcChans = 1:32;

% get wavelets
wavs = is_makeWavelet(foi,lfpFs);

try
    display(['trying...' recName])
    % read in the log file
    fileName = ['/lustre/scr/i/a/iain/logFiles/' recName(1:end-7) '.log'];
    fileID   = fopen(fileName);
    formatSpec = '%f %s %s %f %f %f %f %f %f %f %s %f';
    LogFile = textscan(fileID,formatSpec,'HeaderLines',5,'Delimiter', '\t');
    fclose(fileID);
catch; return;
end

% detect images
cellPic = cellfun(@regexp,LogFile(2),{'Picture'},'UniformOutput',false);
picInd  = cellfun(@isempty,cellPic{:}) == 0;
% detect videos
cellVid = cellfun(@regexp,LogFile(2),{'Video'},'UniformOutput',false);
vidInd  = cellfun(@isempty,cellVid{:}) == 0;

% detect stim onsets
onset  = round((find(diff(triggerData(1,:))==1)/Fs)*lfpFs);
offset = round((find(diff(triggerData(2,:))==1)/Fs)*lfpFs);

% define window
win = round(10*lfpFs);

picOnset  = onset(picInd); picOnset(picOnset < abs(win(1)) | picOnset > size(lfpMat,2)+win(1)) = [];
vidOnset  = onset(vidInd); vidOnset(vidOnset < abs(win(1)) | vidOnset > size(lfpMat,2)+win(1)) = [];
recLength = size(lfpMat,2)/lfpFs;
offset(offset > size(lfpMat,2)-win) = [];

picSamps = [];
vidSamps = [];
graySamps = [];
for ipic = 1:numel(picOnset); picSamps = [picSamps picOnset(ipic):picOnset(ipic)+win]; end
for ivid = 1:numel(vidOnset); vidSamps = [vidSamps vidOnset(ivid):vidOnset(ivid)+win]; end
for iev  = 1:numel(offset); graySamps = [graySamps offset(iev):offset(iev)+win]; end

allChans = [ppcChans pulChans'];
numChans = numel(allChans);
for ichan = 1:numChans
    load([recPath 'spikes\cleanSpk_' num2str(allChans(ichan))])
    spkCell{ichan} = round(index*lfpFs);
end

% initialise matrices
ppcPLV_pic  = nan(numFreqs,numChans);
ppcPLV_vid  = nan(numFreqs,numChans);
ppcPLV_gray = nan(numFreqs,numChans);
pulPLV_pic  = nan(numFreqs,numChans);
pulPLV_vid  = nan(numFreqs,numChans);
pulPLV_gray = nan(numFreqs,numChans);
for f = 1:numel(foi)
    display(['img-video spk-PLV ' recName(1:end-7) ' ' num2str(f) '/' num2str(numFreqs)])
    % PPC first
    C_ppc  = conv2(lfpMat(ppcChans,:),wavs{f},'same'); % convole PPC LFP with wavelet
    angPPC = angle(C_ppc); clear C_ppc
    C_pul  = conv2(lfpMat(pulChans,:),wavs{f},'same'); % convole PPC LFP with wavelet
    angPul = angle(C_pul); clear C_pul
    allAng = vertcat(angPPC,angPul); clear angPPC angPul
    
    nspks = 400;
    nreps = 40;
    for ispk = 1:numChans
        spks     = spkCell{ispk};
        
        % grab angles of spikes for each condition
        picPhase  = allAng(:,picSamps(ismember(spks,picSamps))); picPhase(ispk,:) = nan;
        vidPhase  = allAng(:,vidSamps(ismember(spks,vidSamps))); vidPhase(ispk,:) = nan;
        grayPhase = allAng(:,graySamps(ismember(spks,graySamps))); grayPhase(ispk,:) = nan;
        
        % compute PLV if we have sufficient spikes
        if size(picPhase,2) > nspks
            tmp_pic = nan(nreps,numChans);
            for irep = 1:nreps
                rp_pic  = randperm(size(picPhase,2));
                tmp_pic(irep,:)  = abs(mean(exp(1i*picPhase(:,rp_pic(1:nspks))),2));
            end
            % compute mean of spike-PLV to PPC LFP
            ppcPLV_pic(f,ispk)  = nanmean(mean(tmp_pic(:,1:32)));
            % compute mean of spike-PLV to LP/Pulvinar LFP
            pulPLV_pic(f,ispk)  = nanmean(mean(tmp_pic(:,33:end)));
        end
        % compute PLV if we have sufficient spikes
        if size(vidPhase,2) > nspks
            tmp_vid = nan(nreps,numChans);
            for irep = 1:nreps
                rp_vid  = randperm(size(vidPhase,2));
                tmp_vid(irep,:)  = abs(mean(exp(1i*vidPhase(:,rp_vid(1:nspks))),2));
            end
            % compute mean of spike-PLV to PPC LFP
            ppcPLV_vid(f,ispk)  = nanmean(mean(tmp_vid(:,1:32)));
            % compute mean of spike-PLV to LP/Pulvinar LFP
            pulPLV_vid(f,ispk)  = nanmean(mean(tmp_vid(:,33:end)));
        end
        % compute PLV if we have sufficient spikes
        if size(grayPhase,2) > nspks
            tmp_gray = nan(nreps,numChans);
            for irep = 1:nreps
                rp_gray = randperm(size(grayPhase,2));
                tmp_gray(irep,:) = abs(mean(exp(1i*grayPhase(:,rp_gray(1:nspks))),2));
            end
            % compute mean of spike-PLV to PPC LFP
            ppcPLV_gray(f,ispk) = nanmean(mean(tmp_gray(:,1:32)));
            % compute mean of spike-PLV to LP/Pulvinar LFP
            pulPLV_gray(f,ispk) = nanmean(mean(tmp_gray(:,33:end)));
        end
    end
    save([recPath 'analysis/spkPLV_imagesVideo'],'ppcPLV_pic','ppcPLV_vid','ppcPLV_gray','pulPLV_pic','pulPLV_vid','pulPLV_gray');
end
