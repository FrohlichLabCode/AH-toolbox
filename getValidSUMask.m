function validSUMask = getValidSUMask(animalCode,regionNameSpk,suChnID)
% This function convert SU ID into channel ID and then map onto anatomy and
% exclude units not in proper location
validSUMask = true(size(suChnID));
if ismember(regionNameSpk, {'PFC','PPC','VC'})
    % cortical region all SU are valid
    % no change
elseif ismember(regionNameSpk, {'LPl'})
    % This based on Z:\Individual\Angel\FerretExperiment\ChannelMap.xlsx
    switch animalCode
        case '0171'
            validChnLFP = [3:16];
        case '0179'
            validChnLFP = [3:16]; % same as 0171
        case '0180'
            validChnLFP = [1,3:16]; 
        case '0181'
            validChnLFP = [2:5,12:13]; 
    end
    for iSU = 1:numel(suChnID)
        if ~ismember(suChnID(iSU),validChnLFP) % if not in the right position
            validSUMask(iSU) = false; % set mask to 0
        end
    end
end