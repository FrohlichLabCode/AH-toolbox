
animalCode = '0153';
% addpath('E:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\');
addpath('C:\Users\FrohlichLab\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\');
% pathDir = ['E:\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\' animalCode '\toAnalyze\']; % output save to Dropbox
pathDir = ['C:\Users\FrohlichLab\Dropbox (Frohlich Lab)\Codebase\CodeAngel\Ephys\' animalCode '\toAnalyze\']; % output save to Dropbox
saveRootPath = pathDir;
if ~exist(saveRootPath,'dir'); mkdir(saveRootPath); end
cd(pathDir)

% rawPath = ['E:\' animalCode '\rawData\']; % path to raw data files (J:\ for GPU computer, E:\ for Angel's computer)
rawPath = ['Z:\Ferret Data\0153\electrophysiology\']; % remember to check if the Z drive is connected
files   = dir([rawPath animalCode '_AttentionTask3_11*']); % detect files to load/convert
% Get recording names
cl      = struct2cell(files);
nm      = cl(1,:);
for n = 1:numel(nm); name{n} = nm{n}; end
recNames = unique(name);

% loop through each recording
for irec = 1:1%numel(recNames)
    recName = recNames{irec};
    recPath = [rawPath recName '\'];
    cd(recPath)
    
    if ~exist([pathDir recName],'dir'); mkdir(pathDir,recName); end
    if ~exist([saveRootPath recNames{irec}],'dir'); mkdir(saveRootPath,recNames{irec}); end
    savePath = [saveRootPath recNames{irec} '\'];
    if ~exist([savePath 'lfp'],'dir'); mkdir(savePath,'lfp'); end
    lfpPath = [savePath 'lfp\'];
    if ~exist([savePath 'spikes'],'dir'); mkdir(savePath,'spikes'); end
    spkPath = [savePath 'spikes\'];
    
    %%
    % load in the header infromation for the recording
    is_read_Intan_RHD2000_file(recPath,'info.rhd');
    save([savePath 'frequency_parameters'],'frequency_parameters');
    Fs = frequency_parameters.amplifier_sample_rate; % get sample rate
    %%
    % The following MATLAB code reads a timestamp data file and creates a
    % time vector with units of seconds:
    fileinfo = dir('time.dat');
    num_samples = fileinfo.bytes/4; % int32 = 4 bytes
    fid = fopen('time.dat', 'r');
    t = fread(fid, num_samples, 'int32');
    fclose(fid);
    t = t / frequency_parameters.amplifier_sample_rate; % sample rate from header file
    %%
    % The following MATLAB code reads an amplifier data file and
    % creates an electrode voltage matrix with units of microvolts:
    
    num_channels = length(amplifier_channels); % amplifier channel info from header file
    fileinfo = dir('amplifier.dat');
    num_samples = fileinfo.bytes/(num_channels * 2); % int16 = 2 bytes
    fid = fopen('amplifier.dat', 'r');
    v_mat = fread(fid, [num_channels, num_samples], 'int16');
    fclose(fid);
    v_mat = v_mat * 0.195; % convert to microvolts
    %%
    % The following MATLAB code reads an auxiliary input data file and creates
    % a waveform matrix with units of volts:
    num_channels = length(aux_input_channels); % aux input channel info from header file
    fileinfo = dir('auxiliary.dat');
    num_samples = fileinfo.bytes/(num_channels * 2); % uint16 = 2 bytes
    fid = fopen('auxiliary.dat', 'r');
    v_aux = fread(fid, [num_channels, num_samples], 'uint16');
    fclose(fid);
    v_aux = v_aux * 0.0000374; % convert to volts
    %%
    % The following MATLAB code reads a supply voltage data file and creates
    % a waveform matrix with units of volts:
    num_channels = length(supply_voltage_channels); % supply channel info from header file
    fileinfo = dir('supply.dat');
    num_samples = fileinfo.bytes/(num_channels * 2); % uint16 = 2 bytes
    fid = fopen('supply.dat', 'r');
    v_supply = fread(fid, [num_channels, num_samples], 'uint16');
    fclose(fid);
    v_supply = v_supply * 0.0000748; % convert to volts
    %%
    % The following MATLAB code reads a board ADC input data file and
    % creates a waveform matrix with units of volts:
    try
        num_channels = length(board_adc_channels); % ADC input info from header file
        fileinfo = dir('analogin.dat');
        num_samples = fileinfo.bytes/(num_channels * 2); % uint16 = 2 bytes
        fid = fopen('analogin.dat', 'r');
        adc_data = fread(fid, [num_channels, num_samples], 'uint16');
        fclose(fid);
        adc_data = adc_data * 0.000050354; % convert to volts
        adc_data = adc_data';
        save([savePath 'adc_data'],'adc_data','Fs','-v7.3'); clear adc_data
    catch
        % do nothing
    end
    %%
    %The following MATLAB code reads a board digital input data file and
    % creates vector of 16-bit words:
    fileinfo = dir('digitalin.dat');
    num_samples = fileinfo.bytes/2; % uint16 = 2 bytes
    fid = fopen('digitalin.dat', 'r');
    Din = fread(fid, num_samples, 'uint16');
    fclose(fid);
    
    triggerData = zeros(2, length(Din)); % The number of TTL pulses is hard coded now, but can be soft-coded by looking at the maximum value of Din and figuring out the number of TTL pulses from that maximum value (decomposing the value to a linear combinantion of powers of 2). EN
    ind = find(Din==2); % find time instantences where only port 1 is on
    triggerData(1, ind)=1; % set time instantences where only port 1 is on
    ind = find(Din==4); % find time instantences where only port 2 is on
    triggerData(2, ind)=1; % set time instantences where only port 2 is on
    ind = find(Din==6); % find time instantences where both ports are on
    triggerData(:, ind)=1; % set time instantences where both ports are on
    
    
    
    save([savePath 'triggerData'],'triggerData','Fs','-v7.3'); clear triggerData
    
    %%
    % The following MATLAB code reads a board digital output data file and
    % creates vector of 16-bit words:
    try
        fileinfo = dir('digitalout.dat');
        num_samples = fileinfo.bytes/2; % uint16 = 2 bytes
        fid = fopen('digitalout.dat', 'r');
        digital_out_word = fread(fid, num_samples, 'uint16');
        fclose(fid);
    catch
        % do nothing
    end
    %%
    % define a low pass filter at 300Hz
    downFac = 30;         % factor to downsample raw signal for LFP
    lpFreq  = 300;        % lowpass frequency
    ford    = 4;          % filter order
    lfpFs   = Fs/downFac; % LFP sample rate
    nyqFreq = Fs/2;       % Nyquest frequency of raw signal
    [a,b]   = butter(ford,lpFreq/nyqFreq,'low'); % lowpass filter
    [dr,nr] =  rat(Fs/(Fs/downFac)); % Ratio for downsampling
    lfpMat  = []; % initialize lfp matrix
    %%
    NumChan = size(v_mat, 1);
    % loop through list of channels
    for ichan = 1:NumChan
        ichan
        v = v_mat(ichan, :);
        display(['computing lfp for rec ' recName 'chan ' num2str(ichan)])
        [spkTime,spkWav] = detectSpikes(v,Fs); % detect spikes
        save([spkPath 'spk_' num2str(ichan)],'spkTime','spkWav');
        % compute LFP
        lfp             = filtfilt(a,b,v); % lowpass filter raw signal to obtain LFP
        lfpMat(ichan,:) =  resample(lfp,nr,dr); % downsample data
    end
    save([lfpPath 'lfpMat'],'lfpMat','lfpFs','-v7.3'); clear lfpMat
    %%
    % sometimes we don't have digital data so we need a try/catch to avoid errors
    try
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
end



