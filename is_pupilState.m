function is_pupilState(recPath)
% This function divides the pupil diameter time series into numDiv
% different pupil 'states'. These can then be saved to use for state
% dependent analysis later

% AH modified on 10/21/2018
% changed mean to median for elimThresh to reduce the effect of outlier
% this is used for one eye processing


pup = is_load([recPath 'pupil'],'cpup'); % load pupil data

% divide up pupil data  
numDiv     = 5; % divide into 8 equal (in samples) parts
pupDiv     = pup;
elimThresh = [nanmedian(pupDiv,2)-3*nanstd(pupDiv,0,2) nanmedian(pupDiv,2)+3*nanstd(pupDiv,0,2)]; % eliminate outliers for this analysis
pupDiv(pupDiv<elimThresh(:,1)) = nan; % get rid out stark outliers
pupDiv(pupDiv>elimThresh(:,2)) = nan; % get rid out stark outliers
[counts,centers] = hist(pupDiv,1000); % histogram for pupil diameter
cs = cumsum(counts)/sum(counts); % The cumulative sum for the bin counts for each pupil diameter
for idiv = 1:numDiv
    binMin = (idiv-1)*(1/numDiv); % bin lower limit in percentile of total samples
    binMax = idiv*(1/numDiv); % bin upper limit in percentile of total samples
    binInd = find(cs>=binMin & cs<+ binMax);
    binC   = centers(binInd); % get bin centers
    binLims(idiv,:) = [min(binC) max(binC)]; % define bin limits
    binIndex(idiv,:) = (pupDiv>=min(binC) & pupDiv<=max(binC));
    binCenters{idiv} = binC;
    binCounts{idiv}  = counts(binInd);
end
save([recPath 'pupilStates'],'binLims','binIndex','binCenters','binCounts','-v7.3');
