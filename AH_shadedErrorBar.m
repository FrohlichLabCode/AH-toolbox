function [stats,l1] = AH_shadedErrorBar(xvec, ymat, color, errType)
% l1 is handle of the line
% eg. to display legend: legend([l1 l2],'theta','gamma','Location','east');
% ymat is nChn x nBin (match xvec)
xvec = reshape(xvec,1,[]);
stats.median = nanmedian(ymat,1);
stats.mean   = nanmean(ymat,1);
stats.std    = nanstd(ymat,[],1);
stats.sem    = nanstd(ymat,[],1)/sqrt(size(ymat,1));

if strcmp(errType, 'std') == 1
    dy = stats.std;
elseif strcmp(errType, 'sem') == 1    
    dy = stats.sem;
else
    error('errType has to be ''std'' or ''sem''');
end

h = fill([xvec;flipud(xvec)],[y-yd;flipud(y+dy)],color,'linestyle','none');
set(h,'facealpha',.5); % make translucent
l1 = line(xvec,y,'Color',color);