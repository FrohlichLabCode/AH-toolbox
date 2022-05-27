function [PLI,wPLI,dwPLI] = phaseLagIndex_JR(analyticSignal1,analyticSignal2)
% Code adapted from the FieldTrip software toolbox

% Calculate cross spectrum (which is input to fieldtrip function:
% "ft_connectivity_wpli.m")
% This is a means of calculating the vector which is the phase difference
% between the two signals
crossSpectrum = analyticSignal1 .* conj(analyticSignal2);

% Cross spectrum can also be calculated as:
%   Convert to polar angle and subtract
% angleDiff = angle(analyticSignal1) - angle(analyticSignal2);
%   Convert back into analytic signal
% crossSpectrum = exp(1i*angleDiff);

% Imaginary component of the cross spectrum
imagCrossSpectrum = imag(crossSpectrum);
 
% Phase lag index: ratio of phase difference from 0 to pi vs pi to 2pi
PLI = abs(mean(sign(imagCrossSpectrum)));

% PLI can also be calcualted as:
% binarizePhaseDiff = (imagCrossSpectrum > 0 & imagCrossSpectrum < pi)
% phaseRatio = length(find(binarizePhaseDiff)) / length(binarizePhaseDiff);
% PLI = abs(phaseRatio - 0.5) * 2;
 
% Weighted phase lag index
% sum across imaginary values - essentially larger values have stronger
% connectivity - favors the 90 and 270 degree connections
outsum = sum(imagCrossSpectrum);
% Weight the connections by the absolute value of these values
% Remove the bias towards 90 and 270 degrees by removing magnitude effects
absImagCross = abs(imagCrossSpectrum);
outsumW = sum(absImagCross); % output weight
% Normalize the phase lag index by magnitude of values
%wPLI = outsum ./ outsumW;
wPLI = outsum / outsumW;
% This value can be positive or negative but people take the
% absolute value and ignore the direction of the effect
wPLI = abs(wPLI);
 
% Debiased weighted phase lag index
outssq = sum(imagCrossSpectrum.^2);
%dwPLI = (outsum.^2 - outssq)./(outsumW.^2 - outssq);
dwPLI = (outsum^2 - outssq)/(outsumW^2 - outssq);
% The value will always be negative so take the absolute value
dwPLI = abs(dwPLI);

end