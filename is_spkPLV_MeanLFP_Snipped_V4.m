function [evtSpkPLV,evtSpkAngle,numBins] = is_spkPLV_MeanLFP_Snipped_V4(spkCell,evtTime,analyzeTheseTrials,C_mean,lfpFs,f,regionChn,twin)
% This function computes spike phase locking across the entire recording,
% as well as the time-resolved spike phase locking around saccades. 
% input:  spikeCell - spike times in seconds stored in a cell. Dimensions(numChans/ROIs,evt/trials)
%         sacSamp   - saccades in samples
%         Cmean     - The result of the convolution fo an lfp matrix with a
%                     complex Morlet wavelet. Dimensions(evt/trials,samples)
%         lfpFs     - LFP sample frequency
%         f         - index of frequency
%
% output: sponSpkPLV - the estimated spike phase locking value
%                      computed using a predefined number of randomly 
%                      drawn spikes over various repetitions. Dimensions(numChans,1)
%         sacSpkPLV  - The time-resolved spike-PLV computed around the
%                      timing of saccades. Dimensions(channels,time bins) 
% I.S 2016

display(['processing spike PLV ' num2str(f)])

winSize  = 0.55; % size of window to pool spikes (moving window for time resolved); in seconds
halfWin  = winSize/2;
stepSize = 0.01; % sliding window increment size
nSpks    = 2;  % usually 200 number of spikes to compute spike PLV in each time bin; usually 200
nReps    = 50;  % usually 200; number of times to recompute spike PLV on randomly drawn spikes
totSpkThresh = 10; % usually 1000; number of spikes to consider when computing spontaneous PLV; skip channel if doesn't pass this thresh
binC     = twin(1):stepSize:twin(2); % vector of the center of spike PLV bins
numBins  = numel(binC);
numEvt = numel(evtTime);

lfpWinSamps = (twin(1)*lfpFs:twin(2)*lfpFs);

numAllChans   = size(spkCell,1);
evtSpkPLV  = nan(numAllChans,numBins);
evtSpkAngle = nan(numAllChans,numBins);
evtPhase = cell(numEvt,numBins); % initialise cell to fill with spike phases
skippedChn = [];
for iChn = 1:numAllChans
    numTotSpikes = numel([ spkCell{iChn,analyzeTheseTrials} ]);
    display(['Freq: ' num2str(f) ', Spike PLV chan: ' num2str(iChn) ', Total spks: ' num2str(numTotSpikes)])
    % skip channel if there aren't enough spikes 
    if numTotSpikes < totSpkThresh; skippedChn = [skippedChn iChn]; continue; end
   
    lfpPhase = angle(C_mean);

    %%%%%%% Now compute time-resolved spike PLV in window relative to event
    %%%%%%% time
    
    % find spikes that occured in saccade window
    numSpkPerBin   = nan(numEvt,numBins); % Initialise matrix to count spikes per event and time bin
    
    for iEvt = 1:numEvt
        
        thisEvt = analyzeTheseTrials(iEvt);
        
        spks     = spkCell{iChn,thisEvt}; % in seconds
        spkSamp  = round(spks*lfpFs); % spkSamp are Ca peak samples relative to the Ca TTL
        
        % phase data extracted around window centered at ephys TTL
        lfpPhase_trial = lfpPhase( lfpWinSamps + round( evtTime(iEvt)*lfpFs ) );
        
        for ibin = 1:numBins
           
            
            tmpWin  = [binC(ibin)-halfWin binC(ibin)+halfWin]-twin(1); % time window for saccade and time bin
            % first compare spk times for each bin's time window; if a spk
            % happens in the window, extract the spk sample (lfpFs)
            binSpks = spkSamp(spks>tmpWin(1) & spks<tmpWin(2)); % get spike time in samples that fit in timw window of interest
            % use the spk samples to reach into lfpPhase (lfpFs)
            evtPhase{iEvt,ibin} = lfpPhase_trial(binSpks); % fill cell with spike phases; for each bin, populate vector with spike phases
           
            
            % just note down how many spikes found 
            if isempty(binSpks); numSpkPerBin(iEvt,ibin) = 0; else numSpkPerBin(iEvt,ibin) = numel(binSpks); end % count spikes in time window

        end
    end
    
    % for each time bin, pool phases from each spike. Then compute spike
    % PLV across randomly selected spikes for multiple iterations
    nspksPerBin = sum(numSpkPerBin); % the number of spikes per time bin across all events
    if min(nspksPerBin) < nSpks
        skippedChn = [skippedChn iChn];
        display(['Skipped Chn: ' num2str(iChn) 'Num min spk per bin: ' num2str( min(nspksPerBin)) ]);
        continue % if there aren't enough spikes to compute PLV then skip channel
    else
        display([ 'Num min spk per bin: ' num2str( min(nspksPerBin)) ]);
   
       
        for ibin = 1:numBins
            % display(num2str(ibin))
            % collect all spike phases for the time bin across saccades
            tmpPhase = [];
            for iEvt = 1:numEvt
                tmpPhase = horzcat(tmpPhase,evtPhase{iEvt,ibin});
            end
            % loop through nReps times and compute PLV on nSpks number of
            % randomly drawn spikes
            tmpPLV = nan(1,nReps);
            for irep = 1:nReps
                rp           = randperm(numel(tmpPhase));
                tmpPLV(irep) = (nanmean(exp(1i*(tmpPhase(rp(1:nSpks))))));
            end
            evtSpkPLV(iChn,ibin) = mean(abs(tmpPLV));
            evtSpkAngle(iChn,ibin) = angle(mean(tmpPLV));
        end
    end
   
end

meta.skippedChn = skippedChn;
meta.winSize = winSize;
    

