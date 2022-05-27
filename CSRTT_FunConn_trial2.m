clear
tic

cluster = 1;
skipRec = 1;
linORlog = 2; %freqs of interest: 1=linear 2=log
MedianorPCA = 0; 
animalCode = '0147';


if cluster == 0 %linux use '/', windows matlab can use both '/' and '\'
    addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys'));
    PreprocessDir = ['E:/FerretData/' animalCode '/Preprocessed/'];
    AnalysisDir   = ['E:/FerretData/' animalCode '/Analyzed/'];
    BehavDatDir   = ['E:/FerretData/' animalCode '/behav/'];
    GroupAnalysisDir = ['E:/FerretData/' animalCode '/GroupAnalysis/Spec_missed/'];
elseif cluster == 1
    addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
    PreprocessDir = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Preprocessed/']; % CHANGE FOR KILLDEVIL VS LONGLEAF
    AnalysisDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Analyzed/'];
    BehavDatDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/behav/'];
    GroupAnalysisDir = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/GroupAnalysis/Spec_missed/'];

    %code for initialising parallel computing
    numCore = 24; % USR DEFINE
    myPool = parpool('local',numCore,'SpmdEnabled',false);  
    %Because parfor iterations do not involve interworker communication, disabling SPMD support this way allows the parallel pool to keep evaluating a parfor-loop even if one or more workers aborts during loop execution. 
end

% Different folder name if use median or PCA LFP
if MedianorPCA == 0;     folderSuffix = '_validChns_150f'; %use all valid channels
% elseif MedianorPCA == 1; folderSuffix = '_median'; %use median of all valid channels
% elseif MedianorPCA == 2; folderSuffix = '_PCA';
elseif MedianorPCA == 3; folderSuffix = '_random'; % randomly pick 1 channel
end

fileInfo = dir([PreprocessDir animalCode '_AttentionTask6*']); % detect files to load/convert  '_LateralVideo*'

% loop through each recording
parfor irec = 1:numel(fileInfo)
    try% so that if one record doesn't work, others can still run
    CSRTT_rec_cluster_2(irec,fileInfo,folderSuffix, PreprocessDir, AnalysisDir, BehavDatDir, GroupAnalysisDir,...
        cluster, skipRec, linORlog, MedianorPCA);
    catch
    end
end % end of all records for an animal
delete(myPool)

x = toc;
fprintf('time required =%f sec\n', x);