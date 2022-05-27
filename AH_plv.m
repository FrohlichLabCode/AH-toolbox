function [plv] = AH_plv(xC,yC)
phaseDiff = angle(xC) - angle(yC); % compute phase angle lag between signals (note: not the diff between analytical signals)
if numel(size(xC)) == 2
    % xC, yC are complex matrix of nFOI x nTime
    % Output is nFOI x 1
    plv     = squeeze(abs(nanmean(exp(1i*phaseDiff),2)));
elseif numel(size(xC)) == 3
    % xC, yC are complex matrix of nEvent x nFOI x nTime
    % Output is nFOI x nTime (can then average across time)
    plv     = squeeze(abs(nanmean(exp(1i*phaseDiff),1))); % PLV formula. average across trials (1st dim)
end
end