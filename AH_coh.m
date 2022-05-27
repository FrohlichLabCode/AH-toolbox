function [coh,Cy,icoh] = AH_coh(xC,yC)
% This code will calculate the coherence between two complex signals,
% and output complex coherency Cy as well as amplitude of coherence coh. Can use img
% to get imaginary part of the Cy -- iCoherence.
%
% AH 20210717
if numel(size(xC)) == 2
    % xC, yC are complex matrix of nFOI x nTime
    % Output is nFOI x 1
    Sxy = nanmean(xC.*conj(yC),2); % get coherence of these time columns, average across time, nFOI x 1
    Sxx = nanmean(xC.*conj(xC),2);
    Syy = nanmean(yC.*conj(yC),2);
    Cy  = reshape(Sxy./(sqrt(Sxx.*Syy)),size(xC,1),1); % nFOI x 1, complex number
elseif numel(size(xC)) == 3
    % xC, yC are complex matrix of nEvent x nFOI x nTime
    % Output is nFOI x nTime (can then average across time)
    Sxy = nanmean(xC.*conj(yC),1); % get coherence of these time columns, average across time, nFOI x 1
    Sxx = nanmean(xC.*conj(xC),1);
    Syy = nanmean(yC.*conj(yC),1);
    Cy  = reshape(Sxy./(sqrt(Sxx.*Syy)),size(xC,[2,3])); % nFOI x nTime, complex number
end
coh = abs(Cy); % only get real part, coherence
icoh = imag(Cy); % imaginary part, imaginaryCoherence
end