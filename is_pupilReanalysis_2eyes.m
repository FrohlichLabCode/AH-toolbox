
function is_pupilReanalysis_2eyes(rootPreprocessDir,saveAnalysisDir,condName,excludeTimes,excludeWin)
% modified to accomode 2 eyes sessions (use right eye)
% downsample binIndex from 30000 to 1000Hz to match that of lfpMat 

if ~exist(join(saveAnalysisDir),'dir'); mkdir(join(saveAnalysisDir));end

load([rootPreprocessDir 'pupilStates'])
load([rootPreprocessDir 'lfp/lfpMat']) %change to lfpValid
load([rootPreprocessDir 'saccades'])

if size(pupDiv,1) == 1 
    pup  = pupDiv;
    binIndex = squeeze(binIndex(1,:,:));
    binCounts = binCounts(1,:); % use round bracket to index multiple cells
    binCenters = binCenters(1,:);
    binLims  = binLims(1,:);
    sacSamp = sacSamp{1}; %1 is right eye; based on lfpFs

elseif size(pupDiv,1) == 2 % if there are 2 eyes
    pup = pupDiv(2,:); % use right eye data (better quality)
    %pup      = is_load([rootPreprocessDir 'pupil'],'pupClean'); % load pupil data
    binIndex = squeeze(binIndex(2,:,:)); %2 = right eye data
    binCounts = binCounts(2,:); % use round bracket to index multiple cells
    binCenters = binCenters(2,:);
    binLims  = binLims(2,:);
    sacSamp = sacSamp{1}; %1 is right eye; based on lfpFs
end

condEyeName = [condName '_rPup'];
desacBinIndex = binIndex; 
excludedSamps = excludeTimes*lfpFs;

% remove all time periods around when animals make saccades
remWin = 1*lfpFs;
for isac = 1:numel(sacSamp)
    try
        desacBinIndex(:,(sacSamp(isac) - remWin) : (sacSamp(isac) + remWin)) = 0;
    catch
        % do nothing
    end
end

%numSamps = sum(binIndex,2);%number in 5 bins should be similar since they are 5 equal 20-percentile
numSamps = sum(desacBinIndex,2); 

% bail if there are too few samples per bin to compute func conn
if min(numSamps) < 30*lfpFs; return; end
numDiv   = numel(binCounts);
numChans = size(lfpMat,1);


% Define frequencies of interest. Make sure it's the same as in function
log = 1; % 1 for log spacing of freqs, 0 for linear
numFreqs = 50; % 41; % USR DEFINE
freqBnd = [2 64]; % [1 11]; % USR DEFINE
if log == 1
    foi = logspace(log10(freqBnd(1)),log10(freqBnd(2)),numFreqs); % log spacing
else
    foi = linspace(freqBnd(1),freqBnd(2),numFreqs); % linear spacing
end
wavs         = is_makeWavelet(foi,lfpFs);

pup_powState  = nan(numDiv,numFreqs,numChans);
pup_powNorm   = nan(numDiv,numFreqs,numChans);

desac_powState  = nan(numDiv,numFreqs,numChans);
desac_powNorm   = nan(numDiv,numFreqs,numChans);
desac_PLV       = nan(numDiv,numFreqs,numChans,numChans);
desac_spkPLV    = nan(numDiv,numFreqs,numChans,numChans);

rawFs = 30000;
downFac = 30; 
[dr,nr] =  rat(rawFs/(rawFs/downFac)); % Ratio for downsampling
binIndexDown = logical(resample(double(binIndex)',nr,dr)'); % downsample data, treat column as independent channel
desacBinIndexDown = logical(resample(double(desacBinIndex)',nr,dr)'); 

powSpec = nan(numFreqs,numChans);
for f = 1:numFreqs
    display([rootPreprocessDir 'pupConn ' num2str(f) '/' num2str(numFreqs)])
    % if ~isnan(pupPLV(1,f,1,2)); continue; end
    C = conv2(lfpMat,wavs{f},'same');
    
    % replace excluded epochs with nans
    winRefSamps = excludeWin(1)*lfpFs : excludeWin(2)*lfpFs;
    repMatExcludedSamps = repmat(excludedSamps',[1,numel(winRefSamps)]); % repeat the event samples  
    repMatWinRefSamps = repmat(winRefSamps,[numel(excludedSamps),1]); % repeat the reference window
    allSamps2Rej = round( repMatExcludedSamps + repMatWinRefSamps ); % add the two matrices to get all samples to exclude relative to event samples
    for iChn = 1:size(C,1)
       C(:,allSamps2Rej(:)) = nan;
    end

    % compute mean and standard deviation of power for each channel
    powMean      = nanmean(abs(C).^2,2);
    powStd       = nanstd(abs(C).^2,[],2);
    powSpec(f,:) = powMean;
    
    % Compute the power in different pupil bins
    for ibin = 1:numDiv % loop through previously defined pupil bins
        binData = C(:,binIndexDown(ibin,:)); % get data that corresponds to bin
        powBin                  = nanmean(abs(binData).^2,2);
        pup_powState(ibin,f,:)  = powBin;
        powNorm                 = (powBin-powMean)./powStd;
        pup_powNorm(ibin,f,:)   = powNorm;
        
        % let's do the same after removing saccade epochs
        binData = C(:,desacBinIndexDown(ibin,:)); % get data that corresponds to bin
        powBin                  = nanmean(abs(binData).^2,2);
        desac_powState(ibin,f,:)  = powBin;
        powNorm                 = (powBin-powMean)./powStd;
        desac_powNorm(ibin,f,:)   = powNorm;        
    end
    
    % Now compute other metrics based on pupil derivative
    numReps  = 20;
    numSamps = round(30*lfpFs); % subsample 30 seconds per iteration
    for ibin = 1:numDiv % loop through previously defined pupil bins
        
        % Let's repeat the same analysis again but using the normal pupil
        % bins, with saccade epochs removed
        binData = C(:,desacBinIndexDown(ibin,:)); % get data that corresponds to bin
        
        PLV     = nan(numReps,numChans,numChans);
        % control if we don't have enough data for 30s after removig
        % saccades
        if size(binData,2) > numSamps
            for irep = 1:numReps
                rp      = randperm(size(binData,2)); % random permutation to collect data
                tmpData = binData(:,rp(1:numSamps));
                tmpAng  = angle(tmpData);
                
                for ichan = 1:(numChans-1)
                    for jchan = (ichan+1):numChans
                        ang                   = tmpAng(ichan,:) - tmpAng(jchan,:); % phase difference time series
                        PLV(irep,ichan,jchan) = squeeze(abs(nanmean(exp(1i*(ang))))); % PLV formula
                    end
                end
            end
        end
        % take mean of connectivity measures computed on randomly sampled data
        desac_PLV(ibin,f,:,:)     = squeeze(nanmean(PLV));
        
%         % compute spike phase locking
%         nspks = 200;
%         nreps = 40;
%         for ispk = 1:numChans
%             spks     = is_load([recPath '/spikes/cleanSpk_' num2str(ispk)],'index');
%             pupSpk   = pup(round(spks*lfpFs)); % pupil at time of spik
%             keepInd  = (pupSpk>binLims(ibin,1) & pupSpk<binLims(ibin,2));
%             spkSamp  = round(spks(keepInd)*lfpFs);
%             spkSamp  = spkSamp(ismember(spkSamp,find(desacBinIndex(ibin,:)==1)));
%             for ichan = 1:numChans
%                 spkPhase = angle(C(ichan,spkSamp));
%                 if numel(spkPhase) < nspks; continue; end
%                 tmp_plv = nan(nreps,1);
%                 for irep = 1:nreps
%                     rp = randperm(numel(spkPhase));
%                     tmp_plv(irep) = abs(mean(exp(1i*spkPhase(rp(1:nspks)))));
%                 end
%                 desac_spkPLV(ibin,f,ispk,ichan) = mean(tmp_plv);
%             end
%         end        
    end
    
    
end

save([saveAnalysisDir 'pupFuncConn_' condEyeName],'pup_powState','pup_powNorm',...
        'desac_powState','desac_powNorm','desac_PLV','desac_spkPLV',...
        'foi','excludeTimes','excludeWin','condEyeName','-v7.3');
    
end
