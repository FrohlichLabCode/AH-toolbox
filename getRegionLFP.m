function [regionChn, regionLFP, validChn_1index, newFs] = getRegionLFP(rootPreprocessDir, MedianorPCA, newFs, varargin)
% [regionChn, regionLFP, validChn_1index] = getRegionLFP(rootPreprocessDir,MedianorPCA, newFs)
if nargin == 1
    MedianorPCA = 3; % default using 1 channel
    newFs = [];
elseif nargin == 2
    newFs = [];
end
    
% load lfp
EEG = pop_loadset([rootPreprocessDir 'lfp/lfp_1000fdA.set']);
animalCode = rootPreprocessDir(numel(rootPreprocessDir)-24:numel(rootPreprocessDir)-21);
lfpMat = EEG.data; % resample needs double precision
lfpFs = EEG.srate; %lfpFs  = lfp.Fs; 
lfpChn = EEG.etc; 
% load lfpChn -- OLD
% lfpChn = is_load([rootPreprocessDir 'eeglab_validChn.mat'], 'lfp'); %
% don't use this, some session have more number of channels than EEG
if ~isempty(newFs) % if exist a newFs, downsample
    lfpDownsampled = resample(double(lfpMat)',newFs,lfpFs)';
    lfpInput = lfpDownsampled;
else
    lfpInput = lfpMat;
end
if isempty(newFs); newFs = lfpFs; end % output same lfpFs

numRegion = numel(lfpChn.allChn);
region = getAnimalInfo(animalCode);
regionNames = region.Names;
% correct validChn
lfpChn.clean_channel_mask = logical(lfpChn.clean_channel_mask);
for iRegion = 1:numRegion	
    lfpChn.validChn{iRegion} = lfpChn.allChn{iRegion}(lfpChn.clean_channel_mask((iRegion-1)*16+1:iRegion*16));
end
EEG.etc.validChn = lfpChn.validChn;
EEG.etc.clean_channel_mask = lfpChn.clean_channel_mask;
pop_saveset(EEG,'filepath',[rootPreprocessDir 'lfp/'],'filename',['lfp_1000fdA.set']);


if MedianorPCA == 0 %0=_validChns
    for iRegion = 1:numRegion
        regionChn{iRegion} = lfpChn.validChn{iRegion}; % Pulvinar, PPC, VC
        regionLFP{iRegion} = lfpInput(lfpChn.reorderedChn{iRegion},:); % reordered channel correspond to reordered lfp
        validChn_1index{iRegion} = lfpChn.validChn{iRegion} - 16*(iRegion-1);
    end
elseif MedianorPCA == 4 %4=_validAnaChns (do this since eeglab didn't use anatomical valid channels)
    for iRegion = 1:numRegion
        regionName = regionNames(iRegion); % although only LPl has excluded anatomical channels
        validSUMask = getValidSUMask(animalCode,regionName,lfpChn.validChn{iRegion}-16*(iRegion-1));
        regionChn{iRegion} = lfpChn.validChn{iRegion}(validSUMask); % Pulvinar, PPC, VC
        regionLFP{iRegion} = lfpInput(lfpChn.reorderedChn{iRegion}(validSUMask),:); % reordered channel correspond to reordered lfp
        validChn_1index{iRegion} = lfpChn.validChn{iRegion}(validSUMask) - 16*(iRegion-1);
    end
    
elseif MedianorPCA == 3 %3=_firstChn  
    for iRegion = 1:numRegion
        regionChn{iRegion} = lfpChn.validChn{iRegion}(1);
        regionLFP{iRegion} = lfpInput(lfpChn.reorderedChn{iRegion}(1),:);
        validChn_1index{iRegion} = lfpChn.validChn{iRegion} - 16*(iRegion-1); % for spike-PLV
    end

    if strcmp(animalCode,'0179') % manually pick channels
        optoChns = [7,9,7,0]; % updated on 2/18/2020
        % 0 means unclear which channel is good, use 1st valid channel
        for iRegion = 1:numRegion
            optoChn = optoChns(iRegion)+(iRegion-1)*16;
            newID = lfpChn.reorderedChn{iRegion}(lfpChn.validChn{iRegion}==optoChn);
            if newID>0 % exist this as valid channel then replace
                regionChn{iRegion} = optoChn;
                regionLFP{iRegion} = lfpInput(newID,:);
            end
        end
    elseif strcmp(animalCode,'0180') % manually pick channels
        optoChns = [9,6,4,0]; % updated on 2/18/2020
        for iRegion = 1:numRegion
            optoChn = optoChns(iRegion)+(iRegion-1)*16;
            newID = lfpChn.reorderedChn{iRegion}(lfpChn.validChn{iRegion}==optoChn);
            if newID>0 % exist this as valid channel then replace
                regionChn{iRegion} = optoChn;
                regionLFP{iRegion} = lfpInput(newID,:);
            end
        end
    elseif strcmp(animalCode,'0181') % manually pick channels
        optoChns = [15,13,13,8]; % updated on 2/18/2020
        for iRegion = 1:numRegion
            optoChn = optoChns(iRegion)+(iRegion-1)*16;
            newID = lfpChn.reorderedChn{iRegion}(lfpChn.validChn{iRegion}==optoChn);
            if newID>0 % exist this as valid channel then replace
                regionChn{iRegion} = optoChn;
                regionLFP{iRegion} = lfpInput(newID,:);
            end
        end
    elseif strcmp(animalCode,'0171') % manually pick channels
        optoChns = [4,10,9,0]; % updated[4,13,9,0] on 2/18/2020, [4,10,9,1]on 2/18/2020out of 16 channels per region
        for iRegion = 1:numRegion
            optoChn = optoChns(iRegion)+(iRegion-1)*16;
            newID = lfpChn.reorderedChn{iRegion}(lfpChn.validChn{iRegion}==optoChn);
            if newID>0 % exist this as valid channel then replace
                regionChn{iRegion} = optoChn;
                regionLFP{iRegion} = lfpInput(newID,:);
            end
        end
    end
% if median or PCA, load different lfp, different validChn
elseif MedianorPCA == 1 %1=_median
    if exist([rootPreprocessDir 'lfp/lfpValid.mat']) % older processing pipeline
        lfp = is_load([rootPreprocessDir 'lfp/lfpValid.mat'],'lfp');
        lfpMat = lfp.median;
    else % new pipeline
        for iRegion = 1:numRegion
            lfp.median(iRegion,:) = nanmedian(lfpInput(lfpChn.reorderedChn{iRegion},:),1);
        end
        lfpMat = lfp.median;
        save([rootPreprocessDir 'lfp/lfpValid.mat'],'lfp');
    end
    for iRegion = 1:numRegion %lfp.median is an nChannel by nTimepoint array
        regionChn{iRegion} = iRegion; % Pulvinar, PPC, VC
        regionLFP{iRegion} = lfpMat(iRegion,:); % reordered channel correspond to reordered lfp
        validChn_1index{iRegion} = lfpChn.validChn{iRegion} - 16*(iRegion-1);
    end

elseif MedianorPCA == 2 %2=_PCA;
    EEG = pop_loadset([rootPreprocessDir 'lfp/lfp_1000fdA.set']);
    lfp = is_load([rootPreprocessDir 'lfp/lfpValid.mat'],'lfp');
    lfpMat = lfp.PCA;
    for iRegion = 1:numRegion%lfp.PCA is an nChannel by nTimepoint array
        regionChn{iRegion} = iRegion; % Pulvinar, PPC, VC
        regionLFP{iRegion} = lfpMat(iRegion,:); % reordered channel correspond to reordered lfp
        validChn_1index{iRegion} = lfpChn.validChn{iRegion} - 16*(iRegion-1);
    end
end