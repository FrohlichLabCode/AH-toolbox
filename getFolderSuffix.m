function folderSuffix = getFolderSuffix(MedianorPCA) %0=_validChns; 1=_median; 2=_PCA; 3=_opto1Chn; 4=_validAnaChns;

% Different folder name if use median or PCA LFP
if MedianorPCA == 0;     folderSuffix = '_validChns'; %use all valid channels
elseif MedianorPCA == 1; folderSuffix = '_mdChn'; %use median of all valid channels
elseif MedianorPCA == 2; folderSuffix = '_PCA';
elseif MedianorPCA == 3; folderSuffix = '_opto1Chn'; % randomly pick 1 channel that is optically responsive
elseif MedianorPCA == 4; folderSuffix = '_validAnaChns'; % valid anatomical Chns
end