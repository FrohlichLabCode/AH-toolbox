function [fitresult, gof] = AH_createSineFit(x, y, doPlot)
%CREATEFIT(X,Y)
%  Create a fit.
%
%  Data for 'Sinefit' fit:
%      X Input : x
%      Y Output: y
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 27-Jul-2021 20:21:57
% AH: modified on 7/27/2021 to make initial point close to cosine, add flag
% for plotting

%% Fit: 'Sinefit'.
y = y-median(y); % center y
[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
% Customized eqution % 1 for frequency
ft = fittype( 'a*sin(2*(x-b))+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Normalize = 'on';
opts.StartPoint = [(max(y)-min(y))/2 pi/2 0.1];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

if doPlot
    % Plot fit with data.
    figure( 'Name', 'Sinefit' );
    h = plot( fitresult, xData, yData );
    legend( h, 'y vs. x', 'Sinefit', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % Label axes
    xlabel( 'x', 'Interpreter', 'none' );
    ylabel( 'y', 'Interpreter', 'none' );
    set(gca,'XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})  % to inverse Y axis
    grid on
end

