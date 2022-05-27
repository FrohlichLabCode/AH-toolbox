function [XSmote,YSmote] = AH_SMOTE(XTrain,YTrain,k)
% AH_SMOTE  Synthetic Minority Oversampling Technique. 
% Adapted from a technique to generate synthetic samples as given in: https://www.jair.org/media/953/live-953-2037-jair.pdf
%   Usage:
%   [Xsmote,Ysmote] = AH_SMOTE(XTrain,YTrain,k)
%   
%   Inputs:
%   XTrain: Original dataset features
%   YTrain: Original dataset labels
%   k: number of nearest neighbors (not including itself) to consider while performing
%   augmentation
% 
%   Outputs:
%   XSmote: augmented dataset containing original data as well.
%   YSmote: correponding labels for the augmented dataset
%
%   See also datasample, randsample
% AH created on 2020.5.13

%% plot the bar plot for number of classes
doPlot = 0;
H0counts = histcounts(YTrain); % dont plot, just count
if doPlot == 1
figure
H0 = histogram(YTrain);
ylabel('N of sampels in each class')
xlabel('Class ID')
title('Original imbalance data distirbution')
end
%% number of each classes
labels=YTrain;
classes= categories(YTrain);
for iClass=1:numel(classes)
    classNo(iClass)=numel(find(labels==classes(iClass)));
end

%% Multipliers of sample size for each minority class
[maximumSamples,sampleClass]=max(classNo); % number of maximum samples
for iClass=1:numel(classes)
    % amplify by sample number ratio between classes
    multi(iClass) = floor(maximumSamples / classNo(iClass))-1; % "-1" to subtract the sample itself; use floor so that minority group number doesn't exceed majority group, ceil is fine too, doesn't matter too much
%     % amplification number is how many hundreds (10->100 is 9 times more)
%     % only works if sample size of minority group is around 100;
%     samplediff(ii)=maximumSamples-classNo(ii);
%     multi(ii) = ceil(samplediff(ii)/ 100); 
end

%% Oversample the minority classes, copy the majority class
XSmote=[]; % flattened matrix, nSample x nFeatures
YSmote=[];
for iClass=1:numel(classes)
    classMask = labels==classes(iClass);
    thisN = sum(classMask); % sample size
    thisX = reshape(XTrain(:,:,:,classMask),[],thisN)'; %h x w x 1 x nSample
    thisY = YTrain(classMask);
    XSmoteC = thisX; % each class smote, start with original data
    YSmoteC = AH_cat2num(thisY); % from categorical to double; don't use double, otherwise will change number
    
    if multi(iClass) >= 1 % check if is minority class
        % find k-nearest samples for all points (k+1 include original data
        % point)
        neighborIDs = knnsearch(thisX,thisX,'k',k+1); % row=sample, column=feature, find each X's neighbor in X
        %neighbors is nSample x k matrix  
        for iSample = 1:thisN % For each point, syn multi(ii) new points 
            % retain only N out of k nearest samples, datasample randomly
            % samples allowing for repeats
            neighborID = datasample(neighborIDs(iSample,2:end),multi(iClass),'Replace',true); % find its neighbor index
            thisSeed = thisX(iSample,:);
            thisNeighbor = thisX(neighborID,:);
            % MAIN computation (element-wise)
            thisSyn = (thisNeighbor-thisSeed).*rand(multi(iClass),1) + thisSeed;
            XSmoteC = cat(1, XSmoteC, thisSyn);
            
%             %% Test case: Visualization of an example syn
%             nRand = 3;
%             [foi, tickLoc, tickLabel,~,~] = getFoiLabel(2, 128, 75, 2);
%             xvec = [-4:0.1:0];
%             AH_figure(nRand,3,'SMOTE example')
%             for iRand = 1:nRand
%                 % visualize original seed            
%                 subplot(nRand,3,(iRand-1)*3+1)
%                 imagesc(xvec,1:numel(foi),reshape(thisSeed,size(XTrain,1:2)));
%                 set(gca,'TickDir','out','YTick',tickLoc,'YTickLabel',fliplr(tickLabel))
%                 title(['Label=' num2str(int8(class(ii)))  ' Seed ' num2str(iSample)]);
%                 % visualize this neighbor
%                 subplot(nRand,3,(iRand-1)*3+2)
%                 imagesc(reshape(thisNeighbor(iRand,:),size(XTrain,1:2)));
%                 set(gca,'TickDir','out','YTick',tickLoc,'YTickLabel',fliplr(tickLabel))
%                 title(['Neighbor ' num2str(iRand)]);
%                 % visualize a synthesetic sample
%                 subplot(nRand,3,(iRand-1)*3+3)
%                 imagesc(reshape(thisSyn(iRand,:),size(XTrain,1:2)))
%                 set(gca,'TickDir','out','YTick',tickLoc,'YTickLabel',fliplr(tickLabel))
%                 title(['Synthetic ' num2str(iRand)]);
%             end
        end
        YSmoteC = ones(size(XSmoteC,1),1)*str2num(classes{iClass}); % create Y label
    end    
    XSmote = cat(1,XSmote,XSmoteC); % concatenate along 1st dimension
    YSmote = cat(1,YSmote,YSmoteC);
end
% Convert YSmote to categorical
YSmote = categorical(YSmote);
%% Plot new balanced data histogram
H1counts = histcounts(YSmote);
if doPlot == 1
figure
H1 = histogram(YSmote);
ylabel('N of sampels in each class')
xlabel('Class ID')
title('Balanced data distirbution')
end
%% randomize the data
nSmote = size(XSmote,1);
shuffleID = randperm(nSmote);
XSmote = XSmote(shuffleID,:);
YSmote = YSmote(shuffleID,:);

% resize XSmote from nSample x nFeature
% to h x w x 1 x nSample
XSmote = reshape(XSmote',[size(XTrain,1:3),nSmote]); 
fprintf(['SMOTE k = ' num2str(k) ', nSample = ' num2str(H0counts) ' -> ' num2str(H1counts) '\n']);
end