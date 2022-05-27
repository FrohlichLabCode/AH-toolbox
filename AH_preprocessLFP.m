% This is used after IntanLoad to preprocess LFP data
% Preprocessing steps are:
%       1. reject noise using sub_rejectNoise.m
%       2. select valid channels based on mapping in keepChn.m
%       3. calculate median and PCA for each region
% input: 
%       lfpMat.mat
% output: 
%       lfpDenoised.mat (after step 1)
%       lfpValid.mat (after step 1,2,3)
%
% work for both 5CSRTT and Opto
% AH 6/30/2018

clear
skipRec = 1;

addpath('E:\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\Ephys\');
%animals = {'0181','0180','0171','0173','0172'}; %{'0168','0169'};
animals = {'0187','0188'};
for iAnimal = 1:numel(animals)
    animalCode = animals{iAnimal};
switch animalCode
    case '0172'    
        rootPreprocessDir = ['D:\FerretData\' animalCode '\Preprocessed\'];
    case {'0181','0180','0171','0173'}
        rootPreprocessDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/Preprocessed/'];
    case {'0187','0188'}
        baseDir = ['Z:/Individual/Angel/FerretData/'];
        rootPreprocessDir = [baseDir animalCode '/Preprocessed/'];      
end
cd(rootPreprocessDir)
fileInfo = dir([rootPreprocessDir animalCode '_Opto*']); % detect files to load/convert

% loop through each recording
for irec = 1:numel(fileInfo) %[9:17,19:20,22:28]
    recName = fileInfo(irec).name; %'0168_Opto_010_20180713'
    splitName   = strsplit(recName,'_');
    %if datetime(splitName{4}(1:8), 'InputFormat', 'yyyyMMdd') <= datetime('20180712', 'InputFormat', 'yyyyMMdd'); continue;end

    recPath = [rootPreprocessDir recName '\'];
    if exist([recPath 'lfp\lfpValid.mat'],'file') % check if already preprocessed records
        display([recName ' already processed']);
        if skipRec == 1; continue; end % skip record
    end
    
    display(['processing lfp for rec ' recName]);
    %load lfpMat
    load([recPath 'lfp\lfpMat']);

    % reject noise
    for iChan = 1:size(lfpMat,1)
        lfpDenoised(iChan,:) = sub_rejectNoise(lfpMat(iChan,:),lfpFs,1,1); %(~,~,window to cut out, threshold in uV defined inside function)
        display([recName ' rejected ' num2str(sum(isnan(lfpDenoised(iChan,:)))) ' data points for channel ' num2str(iChan)]);
    end

    % reject channels
    [lfp.validChn, lfp.reorderedChn] = keepChn(recName);
    lfp.reorderedSig = [];

    % calculate median and PCA
    for iRegion = 1:numel(lfp.validChn)
        lfp.reorderedSig = [lfp.reorderedSig; lfpDenoised(lfp.validChn{iRegion},:)]; % signals from valid channels
        lfp.median(iRegion,:) = nanmedian(lfpDenoised(lfp.validChn{iRegion},:),1);
        temp = lfpMat(lfp.validChn{iRegion},:);
        %temp(any(isnan(temp), 2), :) = [];
        [coeff{iRegion},score{iRegion},latent{iRegion},tsquared{iRegion},explained{iRegion},mu{iRegion}] ...
                = pca(temp','VariableWeights','variance', 'Centered',true); % same result
        %PCA can't deal with rows contains NaN; tried 'Row''Complete', doesn't work, return empty value      
        lfp.PCA(iRegion,:) = nanmean(temp.*repmat(coeff{iRegion}(:,1),[1 size(lfpDenoised,2)]),1); % get PCA1 for each region (weighted avg of each channel's lfp, weight=PCA1 coefficient)           

        % scale up lfpPCA
        PCAamp = nanmean(abs(lfp.PCA(iRegion,:))); %VCamp=7
        lfpamp = nanmean(nanmean(abs(temp),2),1);%VCamp=29
        lfp.PCA(iRegion,:) = lfp.PCA(iRegion,:)*lfpamp/PCAamp;%VCamp=29
    end
    lfp.Fs = lfpFs;
    
    display(['saving preprocessed lfp for rec ' recName])
    save([recPath 'lfp\lfpDenoised.mat'], 'lfpDenoised','lfpFs'); %optional
    save([recPath 'lfp\lfpValid.mat'], 'lfp');
    % same channels for spike analysis 
    validChn = lfp.validChn; reorderedChn = lfp.reorderedChn;
    save([recPath 'validChn.mat'], 'validChn', 'reorderedChn'); 
    clear lfp lfpDenoised
end
end