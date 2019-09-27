function matZ = AH_zscore(mat)
% nchn by ntime
[nChn, nt] = size(mat);
mn   = repmat(nanmean(mat, 2), 1, nt);
std  = repmat(nanstd(mat, [], 2), 1, nt);
matZ = (mat - mn)./std;
end