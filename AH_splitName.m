function [animalCode, level, sessionID, date] = AH_splitName(recName)
% works for 5CSRTT file names eg. '0171_Level6b_21_20190419'
splitName  = strsplit(recName,'_');

animalCode = splitName{1};
level      = splitName{2}(6:7);
sessionID  = splitName{3};
date       = splitName{4};


