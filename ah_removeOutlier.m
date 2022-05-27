function [newArray, keptMask] = ah_removeOutlier(array)
% This function will take in an array, calculate its IQR and remove
% outliers that are outlide [Q1-1.5IQR, Q3+1.5IQR]

IQR      = iqr(array);
lowr     = prctile(array,25)-1.5*IQR;
highr    = prctile(array,75)+1.5*IQR;
keptMask = array>lowr & array<highr;
newArray = array(keptMask);
end