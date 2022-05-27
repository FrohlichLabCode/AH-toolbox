function AH_cleanSpikes(recPath)

% recPath = 'E:\Dropbox (Frohlich Lab)\Angel\FerretData\0180\Preprocessed\0180_Level7b_04_20191213\';
% This function is used to clean "fake spikes" that occur in all the
% channels within a defined short time frame (eg. 2 sample)

% AH updated on 2020.1.10, based on is_cleanSpikes.m
% Must satisfy 2 conditions to be count in complimentary Vec
% 1 this spikeTime is a duplicate, 
% 2 this spikeTime has at least 2 more copies within 2 samples (at least 3
% channels have the same spikeTime)

% threshold for defining as co-occurring spikes
maxInterval = 3; % <<-- unit is sample

numChn = numel(dir([recPath 'spikes/spk*.mat']));
load([recPath 'frequency_parameters'])
Fs  = frequency_parameters.amplifier_sample_rate; 
mua = [];
% concatnate all spkTimes into mua
for ichn = 1:numChn
    spkTime = is_load([recPath 'spikes/spk_' num2str(ichn) '.mat'],'spkTime'); 
    mua = [mua round(spkTime*Fs)];
end

[IX,~]  = sort(mua); % sorted spkTimes
diffVec = diff(IX);
cutVec = IX(2:end);
sameVec = cutVec(abs(diffVec) <=maxInterval); % count as same spikes
diffdiffVec = diff(diffVec); % equivalent to minNumChn = 3;
cutVec2 = IX(2:end-1);
dupVec = unique(intersect(cutVec2(abs(diffdiffVec)<=maxInterval),sameVec)); 
%compVec = unique(dupVec);
compVec = unique([dupVec, dupVec-1, dupVec+1]); % also reject spkTime that is 1 sample away

for ichn = 1:numChn
    [spkTime,spkWav] = is_load([recPath 'spikes/spk_' num2str(ichn) '.mat'],'spkTime','spkWav');
    compSpikes   = ismember(round(spkTime*Fs),compVec);
    numRejSpikes = sum(compSpikes); % compute number of rejected spikes
    spkTime(compSpikes==1) = []; % delete co-exist spikes
    rejWav       = spkWav((compSpikes==1),:); % get spkWav for rejected spikes
    normWav      = spkWav((compSpikes==0),:); % get spkWav for accepted spikes
    rejAmp       = min(rejWav,[],2);       
    normAmp      = min(normWav,[],2);
    save([recPath 'spikes/cleanSpk_' num2str(ichn)],'spkTime','numRejSpikes','normWav','rejWav','normAmp','rejAmp');
end
save([recPath 'spikes/corecSpikes'],'compVec');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some debgugging source of corecorded spikes
% 
% numSpike = numel(dir([recPath 'spikes\spk*.mat']));
% mua = [];
% wav = [];
% for ichn = 1:numSpike
%     [index,spkWav] = is_load([recPath 'spikes\spk_' num2str(ichn) '.mat'],'spkTime','spkWav'); 
%     mua = [mua index];
%     wav = [wav ; spkWav];
% end
% 
% [IX,~]  = sort(mua);
% diffVec = diff(diff(IX));
% testVec = IX(2:end-1);
% compVec = unique(testVec(diffVec==0));
% 

% n = [];
% r = [];
% for ichn = 53
%     [spkTime,spkWav]  = is_load([recPath 'spikes\spk_' num2str(ichn) '.mat'],'spkTime','spkWav');
%     compSpikes      = ismember(round(spkTime*Fs),compVec);
%     rejWav          = spkWav((compSpikes==1),:);
%     normWav         = spkWav((compSpikes==0),:);
%     rejAmp          = min(rejWav,[],2);  
%     normAmp         = min(normWav,[],2);
%     r = [r ; rejAmp];
%     n = [n ; normAmp];
% end
% figure();
% subplot(221)
% hist(r,100)
% xlim([-150 0])
% title('corecorded spikes')
% ylabel('count')
% xlabel('amplitude [mV]')
% subplot(223)
% hist(n,100)
% xlim([-150 0])
% title('normal spikes')
% ylabel('count')
% xlabel('amplitude [mV]')
% subplot(222)
% plot(rejWav(1:100,:)'); % plot some sample waveform traces
% hold on
% plot(nanmedian(rejWav,1),'LineWidth',5);
% ylabel('amplitude [mV]');
% xlabel('time [ms]');
% subplot(224)
% plot(normWav(1:100,:)'); % plot some sample waveform traces
% hold on;
% plot(nanmedian(normWav,1),'LineWidth',5);
% ylabel('amplitude [mV]');
% xlabel('time [ms]');