% This code loads in data that have been save in the 'one file per channel'
% INTAN format. A separate function is used to detect spikes (although in 
% the future this should be replaced with the SUA code). 
% I.S. 2015
% AH. 2018: add different animal cases to accomodate for different directory

clear
skipRec = 1;
cluster = 0;
doLFP = 1;
doSpike = 1;
doCleanSpike = 1;
doTrigger = 1;
isACS = 0; %0; whether this is a tACS session = has electrical artifact
addpath(genpath('E:\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\Ephys\'));
animals = {'0201','0187','0188','0186','0182','0185','0158','0183','0180','0181','0179','0171','0173','0172'}; %{'0168','0169'}; % Define which animal we would like to load data from
for iAnimal = 1%:numel(animals)
    animalCode = animals{iAnimal};
    switch animalCode
        case {'0168','0172'} % can add multiple cases
            rootPreprocessDir = ['D:\FerretData\' animalCode '\Preprocessed\']; % output save to Dropbox
            rawPath           = ['D:\FerretData\' animalCode '\rawData\'];
            fileInfo          = dir([rawPath animalCode '_*']); % detect files to load/convert
            CSRTT = 0;
        case {'0173'}
            rootPreprocessDirOld = ['E:\FerretData\' animalCode '\Preprocessed\']; % output save to Dropbox
            rootPreprocessDir = ['E:\FerretData\' animalCode '\Preprocessed_mix\']; 
            rawPath           = ['E:\FerretData\' animalCode '\rawData\'];
            fileInfo          = dir([rawPath animalCode '_Level7*']); % detect files to load/convert
            CSRTT = 1;
        case {'0186','0187','0188'}
            rootPreprocessDir = ['Z:\Individual\Angel\FerretData\' animalCode '\Preprocessed\']; 
            rawPath           = ['Z:\Ferret Data\' animalCode '\ephys\'];
            fileInfo          = dir([rawPath animalCode '_*']); % detect files to load/convert
            CSRTT = 0;
            isACS = 1;    
        case {'0182','0183','0185'}
            rootPreprocessDir = ['Z:\Individual\Angel\FerretData\' animalCode '\Preprocessed\']; 
            rawPath           = ['Z:\Ferret Data\' animalCode '\ephys\'];
            fileInfo          = dir([rawPath animalCode '_*']); % detect files to load/convert
            CSRTT = 0;
        case {'0158'}
            rootPreprocessDir = ['Z:\Individual\Angel\FerretData\' animalCode '\Preprocessed\']; 
            rawPath           = ['Z:\Ferret Data\' animalCode '\ephys\'];
            fileInfo          = dir([rawPath animalCode '_LongLast*']); % detect files to load/convert
            CSRTT = 0;    
        case {'0171','0179','0180','0181'}
            %rootPreprocessDirOld = ['E:\Dropbox (Frohlich Lab)\Angel\FerretData\' animalCode '\Preprocessed\']; % output save to Dropbox
            rootPreprocessDir = ['E:\Dropbox (Frohlich Lab)\Angel\FerretData\' animalCode '\Preprocessed_mix\']; 
            rawPath           = ['Z:\Ferret Data\' animalCode '\ephys\'];           
            if cluster == 1
                rootPreprocessDir = ['/pine/scr/a/n/angelvv/FerretData/' animalCode '/Preprocessed/']; 
                rawPath = rootPreprocessDir;
            end
            fileInfo          = dir([rawPath animalCode '_Level6*']); % detect files to load/convert
            CSRTT = 1;
        case {'0201'}
            rootPreprocessDir = ['E:\Dropbox (Frohlich Lab)\Angel\FerretData\' animalCode '\Preprocessed\']; 
            rawPath           = ['Z:\Ferret Data\' animalCode '\ephys\'];           
            fileInfo          = dir([rawPath animalCode '_Level7*']); % detect files to load/convert
            CSRTT = 1;
    end
    
% loop through each recording
for irec = 2%1:numel(fileInfo)
    recName = fileInfo(irec).name; % Get recording names
    recPath = [rawPath recName '/'];
    if cluster == 0
    if length(recName) < 35 && str2num(recName(1:4))<180 % shorter file name is recorded by block, so process using is_LoadIntanData_Block
        % after animal 0180, file name doesn't have date in it, so only 29 in length, don't want to skip those
        continue
    end

    % detect number of files for each data port
    portAchans = dir([recPath 'amp-A-*']); % PPC
    portBchans = dir([recPath 'amp-B-*']); % Visual cortex 
    portCchans = dir([recPath 'amp-C-*']); % LP/Pulvinar 
    portDchans = dir([recPath 'amp-D-*']); % EEG
    
    % Get the names of all channels and concatenate to get a list of chans
    A = struct2cell(portAchans);
    B = struct2cell(portBchans);
    C = struct2cell(portCchans);
    D = struct2cell(portDchans);
    chanNames = horzcat(A(1,:),B(1,:),C(1,:),D(1,:));
    end
    shortName = recName(1: (length(recName)-14)); % shorten the file name
    if cluster == 1; shortName = recName;end
    savePath = [rootPreprocessDir shortName '/'];
    if ~exist([savePath 'spikes/cleanSpk_1.mat'],'file') || ~exist([savePath 'lfp/lfpMat.mat'],'file')|| ~exist([savePath 'adc_data.mat'],'file')
    %if ~exist([rootPreprocessDir shortName],'dir') || (CSRTT == 1 && ~exist([rootPreprocessDirOld shortName],'dir'))
    elseif skipRec == 1 continue
    end
    
    
    if ~exist([savePath 'lfp/lfpMat.mat'],'file'); AH_mkdir([savePath 'lfp']); doLFP = 1; end % if no lfp then should compute lfp
    lfpPath = [savePath 'lfp/'];
    if ~exist([savePath 'spikes'],'dir'); AH_mkdir([savePath,'spikes']); doSpike = 1; cleanSpike = 1; end
    spkPath = [savePath 'spikes/'];
    
    if (doLFP+doSpike)>0
    % load in the header infromation for the recording
    qf_read_Intan_RHD2000_file(recPath,'info.rhd');
    save([savePath 'frequency_parameters'],'frequency_parameters');
    Fs = frequency_parameters.amplifier_sample_rate; % get sample rate
    
    % define a low pass filter at 300Hz
    downFac = 30;         % factor to downsample raw signal for LFP
    lpFreq  = 300;        % lowpass frequency
    ford    = 4;          % filter order
    lfpFs   = Fs/downFac; % LFP sample rate
    nyqFreq = Fs/2;       % Nyquest frequency of raw signal
    [a,b]   = butter(ford,lpFreq/nyqFreq,'low'); % lowpass filter 
    [dr,nr] =  rat(Fs/(Fs/downFac)); % Ratio for downsampling
    lfpMat  = []; % initialize lfp matrix
    
    % loop through list of channels
    for ichan = 1:numel(chanNames)
        fileName    = chanNames{ichan};
        fileinfo    = dir([recPath fileName]);
        num_samples = fileinfo.bytes/2; %detect file size : int16 = 2 bytes
        fid         = fopen([recPath fileName],'r'); % pointer to file
        v           = fread(fid,num_samples,'int16'); % read data
        fclose(fid); 
        v = v*0.195; % multiply by constant (see INTAN documentation)
        %AH_mkdir([savePath 'raw\']); 
%         display(['saving raw voltage for rec ' recName ' chan ' num2str(ichan)])
%         save([savePath 'raw\rawVol_' num2str(ichan)],'v','Fs');
        % takes 1min per channel only do this for sample session
        
        display(['Computing spikes for rec ' recName ' chan ' num2str(ichan)])
        [spkTime,spkWav] = detectSpikes(v,Fs,isACS); % detect spikes
        save([spkPath 'spk_' num2str(ichan)],'spkTime','spkWav');
        
        % compute LFP
        if doLFP == 1        
        display(['Computing lfp for rec ' recName ' chan ' num2str(ichan)])
        lfp             = filtfilt(a,b,v); % lowpass filter raw signal to obtain LFP
            if ichan > 1
               nSample = size(lfpMat,2);
               tmp = resample(lfp,nr,dr); % downsample data
               lfpMat(ichan,:) = tmp(1:nSample);
            else
               lfpMat(ichan,:) = resample(lfp,nr,dr); % downsample data
            end
        end
    end
    if doLFP == 1
    save([lfpPath 'lfpMat'],'lfpMat','lfpFs','-v7.3'); clear lfpMat
    end
    end
    if doCleanSpike == 1
        display(['Cleaning spikes for rec ' recName])  
        AH_cleanSpikes(savePath); %This function removes co-recorded spikes
    end
    
    % sometimes we don't have digital data so we need a try/catch to avoid errors
    if doTrigger == 1
        Fs = 30000;
    try
        digIn = dir([recPath 'board-DIN*']); % detect digital files
        DI    = struct2cell(digIn);
        DI    = DI(1,:);
        for ichan = 1:numel(DI)
            fileName = DI{ichan};
            fileinfo = dir([recPath fileName]);
            num_samples = fileinfo.bytes/2; % int16 = 2 bytes
            fid = fopen([recPath fileName],'r');
            din = fread(fid,num_samples,'uint16');
            fclose(fid);
            triggerData(ichan,:) =  din; % load digital data into matrix
        end
        save([savePath 'triggerData'],'triggerData','Fs','-v7.3'); clear triggerData
    catch
        % do nothing
    end
    
    % Analogue inputs to INTAN
    try
        ADC = dir([recPath 'board-ADC*']); % detect analogue input files
        ADC_name    = struct2cell(ADC);
        ADC_name    = ADC_name(1,:);
        for ichan = 1:numel(ADC_name)
            fileName = ADC_name{ichan};
            fileinfo = dir([recPath fileName]);
            num_samples = fileinfo.bytes/2; % int16 = 2 bytes
            fid = fopen([recPath fileName],'r');
            v = fread(fid,num_samples,'uint16');
            fclose(fid);
            v = abs(v*0.000050354); % constant obtained from INTAN documentation
            adc_data(ichan,:) = v; % fill ADC matrix
        end
    
        if exist('adc_data')
        save([savePath 'adc_data'],'adc_data','Fs','-v7.3'); clear adc_data
        end
    catch
        % do nothing
    end
    
end
end
end




