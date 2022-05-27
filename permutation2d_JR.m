function [analysisStruct]=permutation2d_JR(data_condBySub,contrastLogic,varargin)

% Optional to include permutation specifications
if ~isempty(varargin)
    permutationOptions = varargin{1};
else
    permutationOptions = struct([]);
end
% Number of iterations to run the permtutation test across
if isfield(permutationOptions,'numIterations')
    numIterations = permutationOptions.numIterations;
else
    numIterations = 1000;
end
% Significance threshold for thresholding across participants
% AND threshold applied to the significance after permutation
if isfield(permutationOptions,'alphaThreshold')
    alphaThreshold = permutationOptions.alphaThreshold;
else
    alphaThreshold = 0.05;
end
% Minimum cluster size to consider for the permutation testing
if isfield(permutationOptions,'minClusterSize')
    minClusterSize = permutationOptions.minClusterSize;
else
    minClusterSize = 30;
end

if isfield(permutationOptions,'logDataFLAG')
    logDataFLAG = permutationOptions.logDataFLAG;
else
    logDataFLAG = 0; % whether take log of data
end

% First two dimensions of the input data are the 2D test dimensions
dim_input = size(data_condBySub);
dim_test = dim_input(1:2);
% Third dimension is the conditions used by "contrastLogic" input
% Fourth dimension is subject to run permutation testing across
numSub = dim_input(3);

% Output is a t-test & difference
analysisOut_t    = NaN(dim_test);
analysisOut_p    = NaN(dim_test);
analysisOut_diff = NaN(dim_test);
analysisOut_mi   = NaN(dim_test);

%% Step 1: cut the data into seprate columns
numColumns = numel(contrastLogic);
colData = zeros([dim_test numSub numColumns]);
% Gather data before running the analysis
% Loop through two columns
for columnIdx = 1:numColumns
    columnLogic = contrastLogic{columnIdx};
    % If there is a null column, then compare to zero
    if columnLogic == 0
        %columnData = zeros([dim_test numSub]); % AH: why is this zero?
        columnData = data_condBySub; % AH 4/23/2021: change to real data
        singleColumn = 1;
    else
        % Average across conditions if multiple
        columnData = squeeze(mean(data_condBySub(:,:,columnLogic,:),3));
        singleColumn = 0;
    end
    % Store this column to run contrast
    colData(:,:,:,columnIdx) = columnData;
end % loop columns
    
% for a single column - run z-transform across data for each participant
if singleColumn
    new_colData = zeros([dim_test numSub]);
    for subIdx = 1:numSub
        subData = squeeze(colData(:,:,subIdx,1));
        subData = reshape(subData,1,numel(subData));
        if logDataFLAG
            subData = log(subData);
        end
        z_subData = zscore(subData);
        new_colData(:,:,subIdx) = reshape(z_subData,dim_test);
    end
    colData = new_colData;
end
    

%% Step 2: loop through each data point and run test
for dim1Idx = 1:dim_test(1)
    for dim2Idx = 1:dim_test(2)
        
        if singleColumn
            col = squeeze(colData(dim1Idx,dim2Idx,:));
            colDiff_all = col;
            colPerc = NaN;
        else
            % Compare two conditions
            col1 = squeeze(colData(dim1Idx,dim2Idx,:,1));
            col2 = squeeze(colData(dim1Idx,dim2Idx,:,2));
            
            if logDataFLAG
                colDiff_all = log(col1) - log(col2);
            else
                colDiff_all = col1 - col2;
            end
            
            % Modulation
            colPerc = mean((col1 - col2) ./ (abs(col1) + abs(col2)));
        end
        % Run t-test
        [~,p,~,tstat]=ttest(colDiff_all);

        % Sometimes column 2 is zeros
        colDiff = mean(colDiff_all);
        
        % Store analysis out
        analysisOut_t(dim1Idx,dim2Idx)    = tstat.tstat;
        analysisOut_p(dim1Idx,dim2Idx)    = p;
        analysisOut_diff(dim1Idx,dim2Idx) = colDiff;
        analysisOut_mi(dim1Idx,dim2Idx)   = colPerc;
    end % loop freq
end % loop time


% Save out tstat, p-values, and decibels for all data channels
analysisStruct_real = struct(...
    't',analysisOut_t,...
    'p',analysisOut_p,...
    'diff',analysisOut_diff,...
    'mi',analysisOut_mi);

% Analysis Structure to save out
analysisStruct = struct(...
    'real',analysisStruct_real);

if numIterations > 0
    % Run permutation testing on data
    iter_maxSize = NaN(1,numIterations);
    iter_maxT = NaN(1,numIterations);
    iter_maxMass = NaN(1,numIterations);
    
    % Loop through iterations
    fprintf('iter ');
    for iterIdx = 1:numIterations
        
        % update the user
        fprintf('%04d',iterIdx);
        
        % Randomize subject assignment
        randSub = randperm(numSub);
        if mod(numSub,2) ~= 0
            extraSub = randi(2) - 1;
        else
            extraSub = 0;
        end
        halfSub = floor(numSub/2) + extraSub;
        randSwapSub = randSub(1:halfSub);
        
        % Gather the t-stat for this iteration
        iterData_t = NaN(dim_test);
        iterData_p = NaN(dim_test);
        
        % Loop through data and run analysis
        for dim1Idx = 1:dim_test(1)
            for dim2Idx = 1:dim_test(2)
                
                if singleColumn
                    col = squeeze(colData(dim1Idx,dim2Idx,:));
                    % Sign flip half of the participants -- this is the key to
                    % create H0 distribution when there is only 1 condition
                    col(randSwapSub) = -col(randSwapSub);
                    colDiff_all = col;
                    tailType = 'right';
                else
                    % Compare two conditions
                    col1 = squeeze(colData(dim1Idx,dim2Idx,:,1));
                    col2 = squeeze(colData(dim1Idx,dim2Idx,:,2));
                
                    % randomly swap half of the subjects
                    tmp = col1;
                    col1(randSwapSub) = col2(randSwapSub);
                    col2(randSwapSub) = tmp(randSwapSub);
                    colDiff_all = col1 - col2;
                    tailType = 'both';
                end
                % run the t-test on random data
                [~,p,~,tstat]=ttest(colDiff_all,zeros(numSub,1),'tail',tailType);
                % Store analysis out
                iterData_t(dim1Idx,dim2Idx) = tstat.tstat;
                iterData_p(dim1Idx,dim2Idx) = p;
                
            end % loop freq
        end % loop time
        
        % Run signifiance test
        sigOptions = struct(...
            'onlyPos',singleColumn);
        [maxSize,maxT,maxMass] = significanceTimeFreq_JR(...
            iterData_t,iterData_p,alphaThreshold,minClusterSize,'test',sigOptions);
        
        % Calculate max size, t, and mass
        iter_maxSize(iterIdx) = maxSize;
        iter_maxT(iterIdx)    = maxT;
        iter_maxMass(iterIdx) = maxMass;
        
        % update the user
        fprintf('\b\b\b\b');
    end % loop iterations
    % Delete the "iter " string preceding the count
    fprintf('\b\b\b\b\b');
    
    % Significance index with alpha = 0.05
    sigIdx = ceil(numIterations * alphaThreshold);
    
    % Sort iterations
    iter_maxSize = sort(iter_maxSize,'descend');
    iter_maxT = sort(iter_maxT,'descend');
    iter_maxMass = sort(iter_maxMass,'descend');
    
    % Get significance threshold
    sig_size = iter_maxSize(sigIdx);
    sig_maxT = iter_maxT(sigIdx);
    sig_mass = iter_maxMass(sigIdx);
    
    % Store output from the permutation testing
    sigStruct = struct(...
        'alpha',alphaThreshold,...
        'minSize',minClusterSize,...
        'size',sig_size,...
        't',sig_maxT,...
        'mass',sig_mass);
    iterStruct = struct(...
        'numIterations',numIterations,...
        'size',iter_maxSize,...
        't',iter_maxT,...
        'mass',iter_maxMass);
    permutationStruct = struct(...
        'sig',sigStruct,...
        'iterations',iterStruct);
    
    % Gather relevant information
    analysisStruct.permutation = permutationStruct;
end % if running permutation analysis
end % end of function