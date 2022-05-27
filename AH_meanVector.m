function [lengths,angles] = AH_meanVector(phaseValue,phaseBins)
% Input: a matrix of n frequency x m phasebins 
% Output: PLV value for each frequency

% 1. Normalize the data for each freq into probability distribution, each row sum to 1
prob = phaseValue./nansum(phaseValue,2);
%Test: allone = nansum(prob,2);

% 2. Calculate mean vector
complexPhase = 1.*exp(1i*phaseBins); % convert phaseBins into unit length vectors
lengths = abs(nansum(prob.*complexPhase,2));
angles = angle(nansum(prob.* complexPhase,2));
end