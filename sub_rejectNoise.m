function sigOut = sub_rejectNoise(sigIn,fs,winLen,rejThresh)
% Note: input sigIn has to be a row or column vector can't be a matrix, because
% of function conv only takes 1D vector
% eg. for fs = 1000Hz, one data point beyond threshold -> 1000 data points
% around become NaN
% use example: xsig = sub_rejectNoise(xser,fs,2,1); 

rejThresh   = 200; % 500uV threshold for rejection -> change to 200uV AH
sigOut      = sigIn;
[b,a]       = butter(4,40/(fs/2),'high'); % define a 4th order highpass filter at 40Hz
hfSig       = filtfilt(b,a,sigIn); % filtyer signal
delWin      = ones(1,round(fs*winLen));                     % window to cut out in sec
delInd      = (abs((hfSig)) > rejThresh);                   % Samples that surpass threshold
delVec      = (conv(double(delInd),delWin,'same') > 0);     % Convolve to smooth window and detect non-zero samples
sigOut(delVec) = nan;                                       % Noisy samples change to NaN's


return