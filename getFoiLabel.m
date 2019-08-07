function [foi, tickLoc, tickLabel,psitickLoc,psitickLabel] = getFoiLabel(lowFreq, highFreq, numFreqs, linORlog)
% [foi, tickLoc, tickLabel,psitickLoc,psitickLabel] = getFoiLabel(2,128,150,2); % lowFreq, highFreq, numFreqs, linORlog)

% Compute ticks for plotting

if linORlog == 1
    foi      = linspace(lowFreq,highFreq,numFreqs); % linear spacing
    psiFreq = foi;
    fois = [lowFreq, 10:10:highFreq];
    tickLabel = string(fois); % generate a string array matches fois {"5","10"...}
    psitickLabel = string([round(psiFreq(1)) fois(2:end-1) round(psiFreq(end))]); % generate a string array for phase slope index
elseif linORlog == 2
    foi   = logspace(log10(lowFreq),log10(highFreq),numFreqs); % log spacing
    psiFreq = foi;
    fois = 2.^(log2(lowFreq):1:log2(highFreq)); %[2 4 8 12 16 32 64 128];
    tickLabel = string(fois);
    psitickLabel = string([round(psiFreq(1)) fois(2:end-1) round(psiFreq(end))]);
end
for fi = 1:numel(fois)
    [bi,bb] = sort(abs(foi-fois(fi)));
    tickLoc(fi) = bb(1);
    [bi,bb] = sort(abs(psiFreq-fois(fi)));
    psitickLoc(fi) = bb(1);
end
return
end