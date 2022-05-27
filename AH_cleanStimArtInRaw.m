% Clean tACS artifact based on template
% AH 2020/10
% recPath = 'Z:\Ferret Data\0187\tmpSpikeSort\0187_arnoldTongue_13_20201022_201022_104746\spikeSort\B\';
% adcPath = 'Z:\Individual\Angel\FerretData\0187\Preprocessed\0187_arnoldTongue_13_20201022\';
% bahavPath = '['Z:/Ferret Data/' animalCode '/behav/tACS/stimFiles/' recName]'
% behavPath = 'Z:\Ferret Data\0187\behav\tACS\stimFiles\0187_arnoldTongue_13_20201022.mat';
function AH_cleanStimArtInRaw(portPath,newPortPath, adcPath, behavPath)
tic
nChn = 16;
stimDur = 90; % sec
fs = 30000; % raw sampling freq

%% Load raw data
fprintf('Loading rawData and cleaning stim artifact: %s \n', portPath(34:end));
rawFile = [portPath 'rawData.dat'];
fid = fopen(rawFile, 'r');
rawData = fread(fid,[nChn,inf],'*int16');
newData = rawData;
fclose(fid);
% figure();plot(tvec, rawData(1,:));xlabel('Time [s]');

%% Load in the tACS waveform data that was fed into the INTAN system
adc_data = is_load([adcPath 'adc_data'],'adc_data'); % adc_data
stim = adc_data(4,:); % 4th channel is stim waveform
tACS_stim = is_load(behavPath,'tACS_stim');
tvec = [1:min(size(rawData,2),numel(stim))]/fs; % sec

stimSmooth = smooth(stim);

% Find stimulation onset by looking for peaks in the first derivative
[pksP,locP] = findpeaks(diff(stim));
[pksN,locN] = findpeaks(-diff(stim));
pkThresh = 0.1; % change from 0.2 to 0.1 to capture 1st stim (AH 2020/10/15)
strLocs = locP(pksP>=pkThresh);
endLocs = locN(pksN>=pkThresh);

% get resting threshold for each channel
% either first 4sec or until the first stim started
restStd = std(single(rawData(:,0.01*fs: min(3.5*fs, strLocs(1)-0.01*fs))),[],2); % need to convert int to single for std
% 
% if numel(strLocs) == numel(endLocs) - 1 % when session stops before a condition ends
%     endLocs(end+1) = numel(stim);
% end
% figure();plot(tvec(2:end),diff(stim));hline(0.1);hline(-0.1);xlabel('Time [s]');
lastsize = 0;
numSamps = stimDur*fs+0.03*fs;
numEv = min(numel(strLocs),numel(tACS_stim));
for iev = 1:numEv % 1s per tACS condition
    fprintf(repmat('\b', 1, lastsize));
    lastsize = fprintf('        stim %d/%d\n', iev, numEv);    
    freq = tACS_stim(iev).freq;
    period = round(1/freq*fs);
    strSamp = strLocs(iev)-0.01*fs; % we want just before stim on to capture the first artifact
    endSamp = min(strSamp+numSamps,numel(stim));
    locs   = strSamp:period:endSamp; % beginning of each cycle
    locsTune = locs;

    evidx   = [strSamp:endSamp];
    stimWav = stim(strSamp:min(strSamp+numSamps,numel(stim))); % get the stimulation waveform, clamp at when recording stops
    stimWav = stimWav-mean(stimWav); % demean
    
    % clean each channel
    for iChn = 1:nChn
        thisData = double(rawData(iChn,evidx)); % only for plotting
        tempWin  = [0,period];
        templates = zeros(numel(locs),period+1);
        % Fine tune cycle alignment based on min value location
        for iloc = 2:numel(locs)-1 % exclude 1st and last sample
            cycleID = locs(iloc):locs(iloc)+period;
            %maxID = find(stimSmooth(cycleID) == max(stimSmooth(cycleID)));
            peakID = find(rawData(iChn,cycleID) == min(rawData(iChn,cycleID)));
            if numel(peakID) > 1; peakID = round(mean(peakID));end % get middle element if not unique
            %if numel(maxID) > 1; maxID = maxID(ceil(end/2));end % get middle element if not unique
            %locsTune(iloc) = cycleID(maxID) - round(1/4*period)-0.01*fs;
            locsTune(iloc) = cycleID(peakID) - round(1/2*period)-0.005*fs;
            
            %{
            figure();plot(stim(cycleID));hold on;plot(stimSmooth(cycleID));legend({'stim','smoothed stim'});
            figure();plot(rawData(iChn,cycleID));
            %}
        end
        locsTune(locsTune+tempWin(2) > size(rawData,2)) = []; % delete locs outside boundary

        % figure();hist(locs-locsTune);title('locs - locsTune');
                
%         % OLD peak-based method
%         cutData = diff(thisData);
%         cutTh = 1.85*restStd(iChn)*max(stimWav)/0.0725; % 0.0725 is the smallest amp for stimWav, baseline
%         cutData(cutData<=cutTh)=0;
%         [pks,locs] = findpeaks(cutData);
%         % Initialization for faster processing
%         tempWin = [-0.03*fs,0.035*fs];
%        templates = zeros(numel(locs)-2,tempWin(2)-tempWin(1)+1);
        
        
        % Collect all templates
        for iloc = 1:numel(locsTune) % don't use first and last to calculate template
            thisLocIdx = [locsTune(iloc)+tempWin(1):locsTune(iloc)+tempWin(2)];
            templates(iloc,:) = double(rawData(iChn,thisLocIdx));
        end
        % Get template for each block
        blockSize = round(1*freq); % cut into 30s blocks
        nBlock = floor(numel(locs)/blockSize); % residuals count towards the last block
        clear mnTemp mdTemp blocks
        blocks = zeros(1,numel(locs)); % all assigned to last block
        for iBlock = 1:nBlock
            if iBlock == 1                
                blocks(((iBlock-1)*blockSize+2) : ((iBlock-1)*blockSize+blockSize)) = iBlock; % 2:30, 31:60 % assign previous blocks
            elseif iBlock < nBlock
                blocks(((iBlock-1)*blockSize) : ((iBlock-1)*blockSize+blockSize)) = iBlock;
            else
                blocks(((iBlock-1)*blockSize+1) : numel(locs)-1) = nBlock;
            end
            mnTemp(iBlock,:) = nanmean(templates(blocks==iBlock,:),1);
            mdTemp(iBlock,:) = nanmedian(templates(blocks==iBlock,:),1);
        end
        finalTemp = mdTemp;
%        
        % modify rawData
        for iloc = 1:numel(locsTune)
            thisLocIdx = [locsTune(iloc)+tempWin(1):locsTune(iloc)+tempWin(2)];
            if iloc == 1 || iloc == numel(locsTune) % set 1st and last cycle to 0             
                newData(iChn,thisLocIdx) = zeros(1,numel(thisLocIdx));
            else
                newData(iChn,thisLocIdx) = rawData(iChn,thisLocIdx) - int16(finalTemp(blocks(iloc),:));
            end
            %{
            figure();subplot(411);plot(rawData(iChn,thisLocIdx))
            subplot(412);plot(templates(iloc,:));
            subplot(414);plot(int16(finalTemp(blocks(iloc),:)))
            subplot(413);plot(newData(iChn,thisLocIdx))
            
            %}
        end

        
        %% Testing block
        %{
        % Plot period-based method
        xLim = [evidx(1)/fs,evidx(0.5*fs)/fs]; % beginning of stim
        xLim = [evidx(end-0.5*fs)/fs,evidx(end)/fs]; % end of stim
        figure();subplot(411);plot(evidx/fs, stimWav);xlim(xLim);ylabel('Stim');
        %vline(locsTune+tempWin(1),'r-'); vline(locsTune+tempWin(2),'k')
        subplot(412);plot(evidx/fs,thisData);xlim(xLim);ylabel('rawData');hold on;
        tempLine = zeros(size(thisData)); tempLine(locsTune(1:end)-strSamp+1) = 100;
        subplot(413);plot(evidx/fs,tempLine);xlim(xLim);
        %vline(locs+tempWin(1),'r-'); vline(locs+tempWin(2),'k')
        
        subplot(413);plot(evidx/fs, diff(thisData));xlim(xLim);xlabel('Time [s]');ylabel('diff(rawData)');hline(cutTh);
        subplot(414);plot(evidx(2:end)/fs, cutData);xlim(xLim);xlabel('Time [s]');ylabel('cutData))');hline(cutTh);        
       
        % Plot templates
        AH_figure(2,3,'');
        subplot(231);plot(templates(1,:)); title('1st cycle');
        subplot(232);plot(templates(2:end-1,:)'); title('middle cycles');      
        subplot(233);plot(templates(end,:)'); title('last cycle'); 
        subplot(234);plot(mnTemp'); title('mean template');
        subplot(235);plot(mdTemp'); title('median template');
        
        % plot after subtraction of template
        AH_figure(3,5,'');
        xLim = [evidx(1)/fs,evidx(0.5*fs)/fs]; % beginning of stim
        subplot(611);plot(evidx/fs, stimWav);xlim(xLim);ylabel('Stim (start)');
        subplot(612);plot(evidx/fs, thisData);xlim(xLim);ylabel('rawData');
        %vline(locsTune+tempWin(1),'r-'); vline(locsTune+tempWin(2),'k');
        subplot(613);plot(evidx/fs, newData(iChn,evidx));xlim(xLim);ylabel('cleanedData');        
        
        xLim = [evidx(end-0.5*fs)/fs,evidx(end)/fs]; % end of stim
        subplot(614);plot(evidx/fs, stimWav);xlim(xLim);ylabel('Stim (end)');
        subplot(615);plot(evidx/fs, thisData);xlim(xLim);ylabel('rawData');        
        %vline(locsTune+tempWin(1),'r-'); vline(locsTune+tempWin(2),'k')
        subplot(616);plot(evidx/fs, newData(iChn,evidx));xlim(xLim);ylabel('cleanedData');
        
        
        
        
        
        
        % Plot threshold-based method
        xLim = [evidx(1)/fs,evidx(0.5*fs)/fs]; % beginning of stim
        xLim = [evidx(end-0.5*fs)/fs,evidx(end)/fs]; % end of stim
        figure();subplot(411);plot(evidx/fs, stimWav);xlim(xLim);ylabel('Stim');
        subplot(412);plot(evidx/fs, thisData);xlim(xLim);ylabel('rawData');
        subplot(413);plot(evidx(2:end)/fs, diff(thisData));xlim(xLim);xlabel('Time [s]');ylabel('diff(rawData)');hline(cutTh);
        subplot(414);plot(evidx(2:end)/fs, cutData);xlim(xLim);xlabel('Time [s]');ylabel('cutData))');hline(cutTh);
        
        % Plot template
        figure();subplot(221);plot(templates');title('allTemp');
        subplot(222);plot(thisData(locs(1)-0.03*fs:locs(1)+0.035*fs));title('egTemp');
        subplot(223);plot(mdTemp);hold on;plot(mdTemp);
        legend({'mnTemp','mdTemp'});
        
        % plot after subtraction of template
        AH_figure(3,5,'');
        xLim = [evidx(1)/fs,evidx(0.5*fs)/fs]; % beginning of stim
        subplot(811);plot(evidx/fs, stimWav);xlim(xLim);ylabel('Stim (start)');
        subplot(812);plot(evidx/fs, thisData);xlim(xLim);ylabel('rawData');
        subplot(813);plot(evidx(2:end)/fs, diff(thisData));xlim(xLim);xlabel('Time [s]');ylabel('diff(rawData)');hline(cutTh);
        subplot(814);plot(evidx/fs, newData(iChn,evidx));xlim(xLim);ylabel('cleanedData');hline(cutTh);
        
        xLim = [evidx(end-0.5*fs)/fs,evidx(end)/fs]; % end of stim
        subplot(815);plot(evidx/fs, stimWav);xlim(xLim);ylabel('Stim (end)');
        subplot(816);plot(evidx/fs, thisData);xlim(xLim);ylabel('rawData');
        subplot(817);plot(evidx(2:end)/fs, diff(newData(iChn,evidx)));xlim(xLim);xlabel('Time [s]');ylabel('cleaned diff');hline(cutTh);
        subplot(818);plot(evidx/fs, newData(iChn,evidx));xlim(xLim);ylabel('cleanedData');
        
        %}
    end
end
rawData = newData;
newFile = [newPortPath 'rawData.dat'];
AH_mkdir(newPortPath);
fid = fopen(newFile,'w');
fwrite(fid,rawData,'int16');
fclose(fid);
end