function tPAC = PAC_time_resolved(lfpMat,fs,event,twin,saveRootPath)
% This function computes a various metrics of functional and effective
% connectivity time-locked to certain events, including:
% 	* Spectral power
% 	* Phase locking value
% 	* Coherence
% 	* Imaginary coherence
% 	* Phase slope index
% 	* Granger causality
% 
% This code also generates a single plot where you can compare all of the
% above-mentioned metrics. 
% Input: - xser  (time series from one brain region)
%        - yser  (time series from another brain region)
%        - fs    (sample rate of time series)
%        - event (vector of event times in seconds)
%        - twin  (time window for spectral analysis, eg twin = [-2 5])
% Output - funcCon (structure containing all information on functional/effective connectivity analyses)
%
% For this script you will need to have the MVGC toolbox in your path. You
% can find the website for the toolbox here: http://users.sussex.ac.uk/~lionelb/MVGC/
% Or you can also just copy the code from the following folder in my
% codebase: C:\Users\FrohlichLab\Dropbox (Frohlich Lab)\Codebase\CodeIain\mainCode\Code\mvgc_v1.0
% I.S. 2017

%% initialize data structure
tPAC = struct;

newFs   = 200;
[nr,dr] =  rat(fs/newFs);
for ichan = 1:size(lfpMat,1) 
    downMat(ichan,:) = resample(lfpMat(ichan,:)',dr,nr)';  % double-checked
end

lowFreqs  = 1:0.25:10;
highFreqs = 15:80;
lowWav    = is_makeWavelet(lowFreqs,newFs);
highWav   = is_makeWavelet(highFreqs,newFs);

for ichan = 1:size(downMat,1)
    testData = downMat(ichan,:); %1 x Timepoints
    lowMat   = nan(numel(lowFreqs),numel(testData));
    highMat  = nan(numel(highFreqs),numel(lowFreqs),numel(testData));

    for lf = 1:numel(lowFreqs)
        lowMat(lf,:) = conv(testData,lowWav{lf},'same'); %nFreq x time
    end

    for hf = 1:numel(highFreqs)
        highTmp = conv(testData,highWav{hf},'same');
        highAmp = abs(highTmp);
        for lf = 1:numel(lowFreqs)
            highMat(hf,lf,:) = conv(highAmp,lowWav{lf},'same');
        end
    end
    
    % Cut out spectral data around events
    tsamps = round(twin*newFs);
    lowmat_event = nan(numel(event),numel(lowFreqs),diff(tsamps)+1); %nTrial x nFreq x nTimepoints
    highmat_event = nan(numel(event),numel(highFreqs),numel(lowFreqs),diff(tsamps)+1);
    for iev = 1:numel(event) % for each trial
        evSamp = round(event(iev)*newFs); % get trial start and end time
        % this may not work if the window looks for samples out of range
        lowmat_event(iev,:,:) = lowMat(:,evSamp+tsamps(1):evSamp+tsamps(2)); % nlf x nhf x nTimepoints
        highmat_event(iev,:,:,:) = highMat(:,:,evSamp+tsamps(1):evSamp+tsamps(2));
    end
    % clear spectral data from memory and compute event-triggered power spectrograms
    clear lowMat highMat

    
    lowAng  = angle(lowmat_event);
    highAng = angle(highmat_event);

    for lf = 1:numel(lowFreqs)
        for hf = 1:numel(highFreqs)
            ang = reshape(lowAng(:,lf,:), size(lowAng,1), size(lowAng,3)) - reshape(highAng(:,hf,lf,:), size(highAng,1),size(highAng,4)); %nTimepoints x 1
            plv(hf,lf,:) = abs(nanmean(exp(1i*ang),1)); %calculate PLV across trials
        end
    end
    tPAC.plv{ichan} = plv;
    
end
tPAC.lowFreqs = lowFreqs;
tPAC.highFreqs = highFreqs;

[filepath,name,ext] = fileparts(saveRootPath);
if ~exist(join(filepath),'dir') 
        mkdir(join(filepath)); end
save(saveRootPath,'tPAC');

return
    
end



