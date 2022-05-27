function AH_getStimOnOffs(recPath,varargin)
%% This script will process adc_data all the way to get stim On Off time
% Outputs stimOns, stimOffs, stim, and correponding plots will be saved into preprocess directory
% AH 2020/11

%{
% Example of use
animalCodes = {'0186','0187','0188','0182','0185'}; % 5 animals all good
for iAnimal = 3%1:numel(animalCodes)
    animalCode = animalCodes{iAnimal};
    rawDir = ['Z:/Ferret Data/' animalCode];
    arnoldFiles = dir([rawDir '/ephys/' animalCode '_arnoldTongue_*']); % raw files
    for irec = 25%1:numel(arnoldFiles)
        recName = arnoldFiles(irec).name(1:end-14);
        recPath = ['Z:/Individual/Angel/FerretData/' animalCode '/Preprocessed/' recName '/'];
        AH_getStimOnOffs(recPath);
    end
end
%}

p=inputParser;
% -- parameters for decomposition -- %
p.addParameter('adcID',4,@isnumeric); % 4th channel is stim
p.addParameter('stimDur',90,@isnumeric); % tACS stim duration is 90s
p.addParameter('lfpFs',1000,@isnumeric); % downsample rate to 1kHz
p.parse(varargin{:});

% Extract varibles from structure fields
% Don't want to save into a function and declare global variable
struct = p.Results;
fieldName = fieldnames(struct);
for ifield=1:length(fieldName)
    eval([fieldName{ifield} '=' ...,
    'struct.(fieldName{ifield});' ]);
end

% Get down-sampled stim
load([recPath 'adc_data'])
nyqFreq = Fs/2;
%lfpFs = Fs/30;
[a,b]   = butter(4,100/nyqFreq,'low');
[dr,nr] =  rat(Fs/lfpFs); % downsample by factor of 30
dat     = filtfilt(a,b,adc_data(adcID,:)); % tACS data are fed into input '4'. 1,2,3 are eye data
stim    = resample(dat,nr,dr);
stim = stim-1.5; % return to baseline
save([recPath 'ACstim.mat'],'stim','lfpFs');
clear adc_data


% Find stimulation onset by looking for peaks in the first derivative
[pks,loc] = findpeaks(diff(stim));
rejInd = (pks<0.2); % change from 0.2 to 0.1 to capture 1st stim (AH 2020/10/15)
pks(rejInd) = [];
loc(rejInd) = [];
stimOns = loc;
% to capture first stim condition
[npks,nloc] = findpeaks(-diff(stim));
stimOffs = nloc(npks>0.2);
stimOffs = stimOffs(stimOffs>0.01*lfpFs);
if stimOffs(1) <= stimOns(1) % missing first stim
    stimOns = [stimOffs(1) - stimDur*lfpFs,  stimOns];
    pks = [pks(1) pks];
end
while stimOffs(end) <= stimOns(end) % for imcomplete last condition
    stimOns(end) = []; % exclude since length will not be the same, so plv calculation will be not fair
    pks(end) = [];
end
if numel(stimOffs) < numel(stimOns) 
%    if stimOffs(end) - stimOns(end) <= (stimDur-10)*lfpFs % more than 10s less than regular duration, then exclude
    stimOns(end) = [];
    pks(end) = [];
%    end
end
% To deal with some weird artifact at end of session
mask = diff(stimOns) < 0.8*stimDur*lfpFs;
mask = [false mask];
stimOns(mask) = [];
pks(mask) = [];
stimOffs(mask) = [];
%&& 
%% Debuging plot
ttl = zeros(1,numel(stim));
tvec = [1:numel(stim)]/lfpFs;
ttl(stimOns) = pks;
ttl(stimOffs) = -pks; % if error, probably stimOn and stimOff are not equal length, need plot for troubleshooting
fig = figure();
recPathPlot = regexprep(recPath(end-29:end-1), '_', ' '); % replace _ with space
subplot(311);plot(tvec,stim);title({recPathPlot; ['stim: ' sprintf('%.1f', numel(stim)/lfpFs/60) ' min']});
subplot(312);plot(tvec(1:end-1),diff(stim),'linewidth',2);title('diff(stim)');
subplot(313);plot(tvec,ttl,'linewidth',2);title({[num2str(numel(stimOns)) ' stimOns'];[num2str(numel(stimOffs)) ' stimOffs']});
xlabel('Time [s]');
save([recPath 'stimOnOffs.mat'],'stim','stimOns','stimOffs','stimDur','pks','lfpFs'); 
savefig(fig, [recPath 'stimOnOffs.fig']);
saveas(fig, [recPath 'stimOnOffs.png']);

end