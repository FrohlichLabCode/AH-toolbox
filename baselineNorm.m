function [normedMat, baseMeanMat,baseStdMat] = baselineNorm(Mat, tvec, baseTwin, doMedian)
% This function takes in a time series or a freq-time matrix and compute
% the baseline normalized matrix. Note that values in the matrix has to be
% real values for median to work. Second dimension of Mat has to be time.

% AH 10/22/2018

% compute baseline mask vector
%tvec = Twin(1):1/fs:Twin(2);
baseMask = tvec>= baseTwin(1) & tvec<= baseTwin(2);
numT = size(Mat,2);

if doMedian == 1
    baseMean = nanmedian(Mat(:,baseMask),2);
else
    baseMean = nanmean(Mat(:,baseMask),2);
end
baseStd = nanstd(Mat(:,baseMask),[],2);

baseMeanMat = repmat(baseMean,1,numT);
baseStdMat  = repmat(baseStd,1,numT);
normedMat   = (Mat-baseMeanMat)./baseStdMat;
