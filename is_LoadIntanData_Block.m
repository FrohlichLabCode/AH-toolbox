% ClusterSwitch = 1; % 1 for analyzing locally
% 

% if ClusterSwitch == 1
%     % paths for local computation
%     addpath(genpath( 'C:\Users\Frohlich Lab\Dropbox (Frohlich Lab)\Codebase\CodeIain')) % '/nas02/home/i/a/iain/Code'))
%     pathDir = ['E:\' animalCode '\rawData\toAnalyze\'];
%     
% elseif ClusterSwitch == 2
%     % % Required for use of KilLDevil and parfor
%     ClusterInfo.setQueueName('week')
%     ClusterInfo.setUserDefinedOptions('-M48')
%     ClusterInfo.setUserDefinedOptions('-R mem128')
%     % % add paths for cluster computation
%     addpath(genpath('/nas02/home/z/h/zhouz/Code/'))
%     pathDir       = [ '/lustre/scr/z/h/zhouz/' animalCode '/'];
%     
%     % code for initialising parallel computing
%     if (matlabpool('size') == 0)
%         cpuInfo = cpuinfo();
%         matlabpool('open',num2str(cpuInfo.NumProcessors));
%     end
%     
% end
% 

% AH updated 2018/9/25

clear
addpath('D:\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\Ephys\');

animals = {'0147'}; %{'0168','0169'}; % Define which animal we would like to load data from
for iAnimal = 1:numel(animals)
    animalCode = animals{iAnimal};


rootPreprocessDir = ['D:\FerretData\' animalCode '\Preprocessed\']; % output save to Dropbox
rawPath           = ['D:\FerretData\' animalCode '\rawData\'];
fileInfo          = dir([rawPath animalCode '*']); % detect files to load/convert

    
if ~exist([rawPath 'tmp'],'dir'); mkdir([rawPath,'tmp']); end
tmpPath = [rawPath 'tmp/'];


for irec = 1:numel(fileInfo)
    recName = fileInfo(irec).name; % Get recording names
    recPath = [rawPath recName '\'];
    
    if length(recName) > 45 % shorter file name is recorded by block, longer file name is recorded by channel
        continue
    end
    if ~exist([rootPreprocessDir recName],'dir'); mkdir([rootPreprocessDir recName]); else continue; end
           
    recPath = [rawPath recName '\'];
    
    recFiles = dir([recPath animalCode '*']);
    for ifile = 1:numel(recFiles)
        is_read_Intan_RHD2000_file(recPath,recFiles(ifile).name);
        save([tmpPath recFiles(ifile).name(1:end-4)],'-v7.3');
        clearvars -except ifile rawPath recFiles fileInfo recName tmpPath irec rootPreprocessDir frequency_parameters animalCode recPath 
    end
    
    matFiles = dir([tmpPath recName '*.mat']);
    data = [];
    adc_data = [];
    triggerData = [];
    for imat = 1:numel(matFiles)
        %[ampData,adc,digData] = is_load([tmpPath matFiles(imat).name],...
          %  'amplifier_data','board_adc_data','board_dig_in_data');
        [ampData,digData] = is_load([tmpPath matFiles(imat).name],...
            'amplifier_data','board_dig_in_data');
        data = horzcat(data,ampData);
        %adc_data = horzcat(adc_data,adc);
        triggerData = horzcat(triggerData,digData);
    end
    
    % grab meta data
    Fs = frequency_parameters.amplifier_sample_rate;
    
    recPath = [rootPreprocessDir recName '\'];
    save([recPath 'frequency_parameters'],'frequency_parameters');
    %save([recPath 'adc_data'],'adc_data','Fs','-v7.3');
    save([recPath 'triggerData'],'triggerData','Fs','-v7.3');
    
    % detect spikes
    if ~exist([recPath 'spikes'],'dir'); mkdir(recPath,'spikes'); else continue; end
    spkPath = [recPath 'spikes\'];
    for ichan = 1:size(data,1)
        display(['computing spikes for rec ' recName 'chan ' num2str(ichan)])
        [index,spikes] = detectSpikes(data(ichan,:),Fs); % index are in ms
        save([spkPath 'spk_' num2str(ichan)],'index','spikes');
    end
    
    % define a low pass filter at 300Hz
    lfpFs = 1000; % CZ switched around to get consistent lfpFs across diff native sampling rates
    downFac = Fs/lfpFs;
    nyqFreq = Fs/2;
    lowPassFreq = 300;
    [a,b] = butter(4,lowPassFreq/nyqFreq,'low');
    [dr,nr] = rat(Fs/(Fs/downFac));
    
    % compute LFP (lfp amplitude is in )
    if ~exist([recPath 'lfp'],'dir'); mkdir(recPath,'lfp'); end
    lfpPath = [recPath 'lfp\'];
    lfpMat = [];
    for ichan = 1:size(data,1)
        display(['computing lfp for rec ' recName 'chan ' num2str(ichan)])
        lfp = filtfilt(a,b,data(ichan,:));
        lfpMat(ichan,:) = resample(lfp,nr,dr);
    end
    
    save([lfpPath 'lfpMat'],'lfpMat','lfpFs');
    
    % delete fiels in temporary folder
    filestodelete = dir([tmpPath '*']);
    for idel = 1:length(filestodelete)
        delete([tmpPath filestodelete(idel).name]);
    end
    clearvars -except ifile rawPath recFiles fileInfo recName tmpPath irec rootPreprocessDir frequency_parameters animalCode recPath 
    
    % end
end
end