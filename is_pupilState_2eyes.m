function is_pupilState_2eyes(recPath)
% This function divides the pupil diameter time series into numDiv
% different pupil 'states'. These can then be saved to use for state
% dependent analysis later

% AH modified on 10/21/2018
% 1. modified to accommodate 2 eyes data, 
% 2. changed mean to median for elimThresh to reduce the effect of outlier 
% 3. added 1st dimension of iEye to output variables -- be careful when use
% those for future analysis


pup = flip(is_load([recPath 'pupil'],'cpup'),1); % load pupil data
eyeNames = {'lPupil', 'rPupil'};

% divide up pupil data  
numDiv     = 5; % divide into 8 equal (in samples) parts
pupDiv     = pup;
elimThresh = [nanmedian(pupDiv,2)-3*nanstd(pupDiv,0,2) nanmedian(pupDiv,2)+3*nanstd(pupDiv,0,2)]; % eliminate outliers for this analysis
pupDiv(pupDiv<elimThresh(:,1)) = nan; % get rid out stark outliers
pupDiv(pupDiv>elimThresh(:,2)) = nan; % get rid out stark outliers

figure();plot(pup');hold on; plot(pupDiv');legend('lPup','rPup','lPupElim','rPupElim');
for iEye = 1:size(pup,1)
    [counts,centers] = hist(pupDiv(iEye,:),1000); % histogram for pupil diameter
    cs = cumsum(counts)/sum(counts); % The cumulative sum for the bin counts for each pupil diameter

    for idiv = 1:numDiv
    binMin = (idiv-1)*(1/numDiv); % bin lower limit in percentile of total samples
    binMax = idiv*(1/numDiv); % bin upper limit in percentile of total samples
    binInd = find(cs>=binMin & cs<= binMax); % find binInd within the percentile eg.0-20%
    binC   = centers(binInd); % get bin centers
    binLims(iEye,idiv,:) = [min(binC) max(binC)]; % define bin limits
    binIndex(iEye,idiv,:) = (pupDiv(iEye,:)>=min(binC) & pupDiv(iEye,:)<=max(binC));
    binCenters{iEye,idiv} = binC;
    binCounts{iEye,idiv}  = counts(binInd);
    end
end
save([recPath 'pupilStates'],'pupDiv','eyeNames', 'binLims','binIndex','binCenters','binCounts','-v7.3');
