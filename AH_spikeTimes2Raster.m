% This script turns a table of spike times into 0,1 encoding raster plot
% input table is a matrix of cells, n unit by n trials
% spike times are already zeroed to aligned with beginning of time window
% AH 2020/2/13

function raster = AH_spikeTimes2Raster(table, twin, timeReso)

[nU, nTrial] = size(table);
for iU = 1:nU
    raster{iU,1} = zeros(nTrial, [twin(2)-twin(1)]*timeReso);
    for iTrial = 1:nTrial
        ID = round(table{iU,iTrial}*timeReso);
        ID(ID <1) = 1; % get rid of 1 index
        ID(ID >[twin(2)-twin(1)]*timeReso) = [twin(2)-twin(1)]*timeReso; % clamp value exceed the window
        raster{iU,1}(iTrial,ID) = 1; % turn spike time sample into 1
    end
end
end
