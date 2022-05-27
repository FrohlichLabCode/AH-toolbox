function [spkTime,spkWav] = detectSpikes(rawData,Fs,stdmin,tACS)
% This function detects negatively deflecting spikes in broadband data. 
% Inputs:  rawData - broadband (or hp filtered) signal recorded invasively
%          Fs      - Sample rate of input signal
% Outputs: spkTime - Vector of the timestamp of each detected spike
%          spkWav  - Matrix of spike waveforms
% I.S. 2016
% A.H. 20180714 added outlier_clamp, adjusted stdmin from -4 to -3.5
% A.H. 20200114 changed outlier_clamp to 2std(rawData), change to median
% threshold instead of std threshold (Quian Quiroga,2004)
if nargin < 4
    tACS = 0;
end

w_pre         = round(0.001*Fs);  % 1ms number of pre-event data points stored (def. 20)
w_post        = round(0.0025*Fs); % 2.5ms number of post-event data points stored (def. 44)
stdmin        = -3;               % minimum threshold (def. 5), change to -3.5 by visual inspection, AH 20190112, change to 3 AH 20200114
detect_fmin   = 300;              % high pass filter for detection 
detect_fmax   = 5000;             % low pass filter for detection 
min_ref_per   = round(0.0015*Fs); % detector dead time (in ms)

% define filter
ford  = 4; % define filter order - this may have to be higher if there are large amplitude low-frequency fluctuations in the data
[b,a] = butter(ford,[detect_fmin detect_fmax]/(Fs/2),'bandpass');

% filter data
hpFiltered = filter(b,a,rawData);
if tACS == 1
outlier_clamp = 4.5*std(hpFiltered(0.01*Fs:4*Fs)); % to avoid outlier deviate std calculation
restmax = nanmax(abs(hpFiltered(0.01*Fs:4*Fs)));
outlier_clamp = min(restmax, outlier_clamp);
outlier_ID = find(abs(hpFiltered)>=outlier_clamp);
k = 0.0005*Fs*5; % AH: from 2020/10/21 note
outlier_ID = AH_expandByK2(outlier_ID,k,[1,numel(hpFiltered)]);
hpFiltered(outlier_ID) = 0; %~0 replace outliers with 0
end
% define threshold (2 methods)
%thresh = stdmin * std(hpFiltered); % 1. standard deviation
thresh = stdmin * median(abs(hpFiltered))/0.6745; % 2. median absolute deviation (Quian Quiroga,2004) AH 1/13/2020

%{
%% Troubleshoot plot
tMask = 1:Fs*60*10;
tvec = [1:tMask(end)]/Fs;
fig = AH_figure(2,4,'');
subplot(2,1,1)
plot(tvec,rawData(tMask));
title(['0187 arnoldTongue 04 20201008 Ch1']);
xlabel('Time [s]'); ylabel('Raw signal');
%xlim([108.6,108.8]);
subplot(2,1,2)
plot(tvec,hpFiltered(tMask)); 
hline(outlier_clamp)
hline(-outlier_clamp)
%title(['0187 arnoldTongue 04 20201008 Ch1']);
xlabel('Time [s]'); ylabel('HpFiltered (300-5000Hz) signal');

% Full tvec
tMask = 1:Fs*60*10;
tvec = [1:numel(hpFiltered)]/Fs; % too long, can't save
fig = AH_figure(2,4,'');
subplot(2,1,1)
plot(tvec,rawData);
title(['0187 arnoldTongue 04 20201008 Ch1']);
xlabel('Time [s]'); ylabel('Raw signal');
%xlim([5338.43,5338.45]);
subplot(2,1,2)
plot(tvec,hpFiltered); 
%title(['0187 arnoldTongue 04 20201008 Ch1']);
xlabel('Time [s]'); ylabel('HpFiltered (300-5000Hz) signal');
hline(thresh);
%xlim([5338.1,5338.5]);
%hline(outlier_clamp); % for AT sessions, the clamp is very far away from
%data, not really exclude anything
% saveas(fig, ['Z:\Individual\Angel\FerretData\0187\sampleHP_resting.fig']);
% saveas(fig, ['Z:\Individual\Angel\FerretData\0187\sampleHP_resting.png']);
savefig(fig, ['Z:\Individual\Angel\FerretData\0187\sampleHP_arnoldTongue.fig'],'compact');
saveas(fig, ['Z:\Individual\Angel\FerretData\0187\sampleHP_arnoldTongue.png']);
%}

%% extract spike indices
spikeInds(1,:) = find( hpFiltered(1:end-1) > thresh & hpFiltered(2:end) <= thresh);

% delete spikes that occur during dead-time
spikeInds(find(diff(spikeInds) < min_ref_per) + 1) = [];
spikeInds(spikeInds<(3*Fs)) = [];
spikeInds(spikeInds>(numel(rawData)/Fs-0.1)*Fs) = [];

spkTime = spikeInds/Fs; % convert sample to time

% extract spike waveforms
spkMat  = repmat(spikeInds',1,numel(-w_pre:w_post));
subMat  = repmat((-w_pre:w_post),numel(spikeInds),1);
spkInd  = spkMat + subMat;
spkSamp = reshape(spkInd',[1 numel(spikeInds)*numel(-w_pre:w_post)]);

spkWav  = reshape(hpFiltered(spkSamp),[numel(-w_pre:w_post) numel(spikeInds)])';
%spkWav  = reshape(rawData(spkSamp),[numel(-w_pre:w_post) numel(spikeInds)])';
% spike waveform from rawData looks noisy, hpfiltered is better --
% see note AH 1/13/2020

%% plot thresholded data
% figure
% tvec = [1:1/Fs:length(hpFiltered)/Fs];
% subplot(211)
% plot(hpFiltered)
% hold on
% tmp = ones(1,length(hpFiltered)).*thresh;
% title('HP filtered signal thresholded');
% plot(tmp);xlabel('Sample@30kHz');ylabel('Amplitude');
% subplot(212)
% plot(hpFiltered)
% hold on
% tmp = ones(1,length(hpFiltered)).*thresh;
% plot(tmp);xlim([1e7,1.05*1e7]); 
% xlabel('Sample@30kHz');ylabel('Amplitude');
% 
% % plot spike waveform
% AH_figure(3,3,'waveform');
% for iSpk = 1:100  
%     subplot(10,10,iSpk)
%     plot(spkWav(iSpk,:)); title(['spk# ' num2str(iSpk)]);
% end


