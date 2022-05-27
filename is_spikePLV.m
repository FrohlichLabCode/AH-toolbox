function phaseSpike = is_spikePLV(freq,recPath,varargin)
% This code takes spike times and computes the spike phase locking value (PLV)
% from complex fourier spectra data computed using ft_freqanalysis.m
% I.S. 2012
% inputs: freq     - fieldtrip structure delivered from ft_freqanalysis
%                    function
%         cfg      - fieltrip config file
%         spikes   - structure size 1 x no. units containing spike times in
%                    seconds
%         timeDiff - time difference between spike and LFP recording onset
%
% outputs: phaseSpike - structure containing PLVs bla bla

% do not calculate PLV for ECOG electrodes
if regexp(freq.label{:},'EU64L'); phaseSpike = []; return; end
    
if isempty(varargin)
    % all conditions the same
    ExpStruct     = [];
    numConditions = 1;
    trialIndex    = 1:size(freq.fourierspctrm,1);
else
    % define trial indexes for all stim conditions
    ExpStruct     = varargin{1};
    numConditions = size(ExpStruct.TrialMtx.mat,1);
    trialCell     = struct2cell(ExpStruct.trials);
    trialCell     = trialCell(3,:);
    for cond = 1:numConditions
        trialIndex(cond,:) = find(cell2mat(cellfun(@(x) x == cond,trialCell,'UniformOutput',0)) == 1);
    end
end

% if spikes are empty then load them in
if length(varargin) == 1
    spkName = ['spikes' freq.label{:}(4:end)];
    load([recPath 'spikes/' spkName '.mat'],'index');
    spikes.times = index/1000;
else
    spikes.times = varargin{2}/1000;
end

trl               = freq.cfg.trl;
trlTime           = freq.time; % 0ms = trial onset
nTrials           = size(freq.cfg.trl,1);
winSize           = 0.1; % sliding window for PLV analysis (ms)
stepSize          = unique(diff(round(trlTime*1000)))/1000;
halfWidth         = winSize/stepSize/2;
numIterations     = 1000;
SRate             = freq.fsample;

for icond = 1:numConditions
    % initialise cell structure to store spike phases
    phaseSpike(icond).rawPhase = cell(1,length(trlTime));
    spkTimes = spikes.times; % get spike times
    
    % loop through trials and place phase angles for each spike into
    % timebins
    for itrial = trialIndex(icond,:)
        tlims  = trl(itrial,1:2)/SRate;
        evTime = tlims(1)-trl(1,3)/SRate;
        binVec = evTime + trlTime;
        spkInd = find(spkTimes>tlims(1) & spkTimes<tlims(2)); % find spikes in trial
        for ispk = spkInd % loop through spikes
            spkT       = spkTimes(ispk); % spike time
            spkmod     = abs(binVec-spkT);
            [B,IX]     = sort(spkmod); % find timebin closest to spike
            spkbin     = IX(1);
            spkang     = squeeze(angle(freq.fourierspctrm(itrial,1,:,spkbin))); % get ang
            phaseSpike(icond).rawPhase{spkbin}(:,end+1) = spkang;
        end
        testing(itrial) = length(spkInd);
    end
    % detect number of spikes in each timebin
    for k = 1:size(phaseSpike(icond).rawPhase,2)
        numspks(k) = size(phaseSpike(icond).rawPhase{k},2);
    end
    numconv = conv(numspks,ones(1,winSize/stepSize),'same');
    numconv(1:winSize/stepSize)       = nan;
    numconv(end-winSize/stepSize:end) = nan;
    maxspks = nanmin(numconv); % number of spikes to use for defined winSize
    
    % intitialise PLV and angle matrices
    PLVmat   = nan(length(freq.freq),length(numconv));
    angMat   = nan(length(freq.freq),length(numconv));
    startbin = winSize/stepSize;
    stopbin  = length(numconv)-winSize/stepSize;
    
    % loop through time bins and measure PLV
    for ibin = startbin:stopbin
        % concatinate spike phases from adjacent timebins
        pMat = horzcat(phaseSpike(icond).rawPhase{ibin-halfWidth:ibin+halfWidth});
        for rep = 1:numIterations % repeat computation and take average
            randInd            = ceil(rand(1,maxspks)*size(pMat,2));
            inpMat             = pMat(:,randInd);
            [X Y]              = pol2cart(inpMat,ones(size(inpMat)));
            Xm                 = nanmean(X,2);
            Ym                 = nanmean(Y,2);
            [ang L]            = cart2pol(Xm,Ym);
            PLVmat(:,ibin,rep) = L;
            angMat(:,ibin,rep) = ang;
        end
    end
    phaseSpike(icond).PLV     = nanmean(PLVmat,3);
    phaseSpike(icond).ang     = nanmean(angMat,3);
    phaseSpike(icond).maxspks = maxspks;
    % save stimulus parameters
    if ~isempty(ExpStruct)
        phaseSpike(icond).evInfo.Auditory = ExpStruct.Stimuli.Auditory;
        phaseSpike(icond).evInfo.Visual   = ExpStruct.Stimuli.Visual;
        phaseSpike(icond).evInfo.AVcommon = ExpStruct.Stimuli.AVcommon;
        phaseSpike(icond).evInfo.params   = ExpStruct.TrialMtx.params;
        phaseSpike(icond).evInfo.mat      = ExpStruct.TrialMtx.mat(icond,:);
    end
end
