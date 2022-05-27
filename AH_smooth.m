function [X2,Y2,smoothedData] = AH_smooth(data, interval)
%// Define integer grid of coordinates for the above data
[X,Y] = meshgrid(1:size(data,2), 1:size(data,1));

%// Define a finer grid of points
[X2,Y2] = meshgrid(1:interval:size(data,2), 1:interval:size(data,1));

%// Interpolate the data and show the output
smoothedData = interp2(X, Y, data, X2, Y2, 'linear');

