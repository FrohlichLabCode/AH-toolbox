function [sfc,sta,freq] = spikeFieldCoherence(spikes,lfp,Fs,win)

% inputs:   spikes - spike times in seconds
%           lfp    - lfp time series
%           Fs     - sample frequency
%           win    - time window to calculate SFC in seconds
%
% outputs   sfc    - spike field coherence
%           sta    - spike triggered average lfp
%           freq   - frequency resolution of spike field coherence    

% We first have to eliminate any spikes that might cause windowing errors
recLength                      = length(lfp)/Fs;
spikes(spikes <= win)           = [];
spikes(spikes >= recLength-win) = [];

% compute number of spikes 
nspk = numel(spikes);

% convert spike times to LFP samples
spikeSamp = round(spikes*Fs);

% compute STA LFP
win    = [-round(win*Fs) round(win*Fs)];
lfpMat = nan(nspk,diff(win)+1);
for ispk = 1:nspk
    lfpMat(ispk,:) = lfp(spikeSamp(ispk)+win(1):spikeSamp(ispk)+win(2));
end

% compute spike triggered average LFP
sta = mean(lfpMat);

% compute the fft of the spike triggered average LFP
L      = length(sta);
n      = 2^nextpow2(L);
Y      = fft(sta,n);
Y      = Y(:,1:n/2+1); % one sided
staPow = Y.*conj(Y)/n; % compute power 

% compute the fft of the individual spike triggered lfp segments
lfpFFT = fft(lfpMat,n,2);
lfpFFT = lfpFFT(:,1:n/2+1); % one sided
lfpPow = lfpFFT.*conj(lfpFFT)/n; % compute power

% compute spike field coherence by normalising STA power with mean power
sfc = staPow./mean(lfpPow);

% compute the frequency 
freq   = (0:n-1)*(Fs/n);
freq   = freq(:,1:n/2+1); % one sided

return