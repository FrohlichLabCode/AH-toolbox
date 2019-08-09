function [regionChn, regionLFP, validChn_1index] = getRegionLFP(rootPreprocessDir, MedianorPCA, newFs, varargin)
if nargin == 1
    MedianorPCA = 3; % default using 1 channel
    newFs = [];
elseif nargin == 2
    newFs = [];
end
    
% load lfp
EEG = pop_loadset([rootPreprocessDir 'lfp/lfp_1000fdA.set']);
lfpMat = EEG.data; % resample needs double precision
lfpFs = EEG.srate; %lfpFs  = lfp.Fs; 
if ~isempty(newFs) % if exist a newFs, downsample
    lfpDownsampled = resample(double(lfpMat)',newFs,lfpFs)';
    lfpInput = lfpDownsampled;
else
    lfpInput = lfpMat;
end

% load lfpChn
lfpChn = is_load([rootPreprocessDir 'eeglab_validChn.mat'], 'lfp');
numRegion = numel(lfpChn.allChn);

if MedianorPCA == 0 %0=_validChns
    for iRegion = 1:numRegion
        regionChn{iRegion} = lfpChn.validChn{iRegion}; % Pulvinar, PPC, VC
        regionLFP{iRegion} = lfpInput(lfpChn.reorderedChn{iRegion},:); % reordered channel correspond to reordered lfp
        validChn_1index{iRegion} = lfpChn.validChn{iRegion} - 16*(iRegion-1);
    end
    
elseif MedianorPCA == 3 %3=_firstChn
    for iRegion = 1:numRegion
        regionChn{iRegion} = lfpChn.validChn{iRegion}(1);
        regionLFP{iRegion} = lfpInput(lfpChn.reorderedChn{iRegion}(1),:);
        validChn_1index{iRegion} = lfpChn.validChn{iRegion} - 16*(iRegion-1);
    end

% if median or PCA, load different lfp, different validChn
elseif MedianorPCA == 1 %1=_median
    lfp = is_load([rootPreprocessDir 'lfp/lfpValid.mat']);
    lfpMat = lfp.median;
    for iRegion = 1:numRegion %lfp.median is an nChannel by nTimepoint array
        regionChn{iRegion} = iRegion; % Pulvinar, PPC, VC
        regionLFP{iRegion} = lfpMat(iRegion,:); % reordered channel correspond to reordered lfp
        validChn_1index{iRegion} = lfp.validChn{iRegion} - 16*(iRegion-1);
    end

elseif MedianorPCA == 2 %2=_PCA;
    lfp = is_load([rootPreprocessDir 'lfp/lfpValid.mat']);
    lfpMat = lfp.PCA;
    for iRegion = 1:numRegion%lfp.PCA is an nChannel by nTimepoint array
        regionChn{iRegion} = iRegion; % Pulvinar, PPC, VC
        regionLFP{iRegion} = lfpMat(iRegion,:); % reordered channel correspond to reordered lfp
        validChn_1index{iRegion} = lfp.validChn{iRegion} - 16*(iRegion-1);
    end
end