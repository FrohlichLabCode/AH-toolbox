function [yscale] = AH_sineFitting(data,x)
% data is nRow by 1 period of data
% Adapted from https://www.mathworks.com/matlabcentral/answers/482647-how-can-i-fit-data-to-a-sine-curve
% and https://www.mathworks.com/matlabcentral/answers/121579-curve-fitting-to-a-sinusoidal-function
% This code will find the best sine wave to fit the data for each row, then extract the amplitude of
% this sine wave.
% 
% for iRow = 1:size(data,1)
%     y = data(iRow,:);
%     % Get estimates for initial fitting
%     yu = max(y);
%     yl = min(y);
%     yr = (yu-yl);                               % Range of ‘y’
%     yz = y-yu+(yr/2);                           % Center y around 0
%     zx = x(yz .* circshift(yz,[0 1]) <= 0);     % Find zero-crossings
%     per = 2*mean(diff(zx));                     % Estimate period
%     ym = mean(y);                               % Estimate offset
%     % Create a model
%     fit = @(b,x)  b(1).*(sin(2*pi*(x - b(2))))+b(3);           % Function to fit (b(2) is shift, b(3) is scale, b(1) is yscale
%     fcn = @(b) sum((fit(b,x) - y).^2);                              % Least-Squares cost function
%     s = fminsearch(fcn, [yr;  per;  -1;  ym]);                      % Minimise Least-Squares
%     xp = linspace(min(x),max(x));
%     yscale(iRow,:) = s(1);
%     
%     % visualize solution
%     figure(1)
%     plot(x,y,'b',  xp,fit(s,xp), 'r')
%     set(gca,'XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})  % to inverse Y axis
% 
%     grid
% end

% Use curvefitting toolbox
% https://www.mathworks.com/products/curvefitting.html
a = cftool(x,data(iRow,:));

% ft = fittype('sin((x - shift)/xscale)*yscale','coefficients',{'shift','xscale','yscale'});
% ft(shift,xscale,yscale,x) = sin((x - shift)/xscale)*yscale;
%      
% % Fit a model
% mdl = fit(X,Y,ft,'startpoint',[shiftguess,xscaleguess,yscaleguess]);
end