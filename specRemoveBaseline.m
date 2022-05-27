function specRemovedBaseline = specRemoveBaseline(tvec, foi, spec, twinBaseline)
nTrials = size(spec,1);
twinBaseline = [-2,-1]; % baseline time window
% normalise to pre saccade firing rate
preBins = (tvec>twinBaseline(2) && tvec<twinBaseline(2)); % baseline bins
baselineMean  = mean(spec(preBins,:,:)),); % mean on the time axis
baselineSTD   = std(spec(preBins,:,:));    % std on the time axis
specRemovedBaseline = (spec-baselineMean)/baselineSTD; % Spike z score
