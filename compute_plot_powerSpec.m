function [pow, foi, tickLoc, tickLabel] = compute_plot_powerSpec(mat,figHandle,numRow, numCol, subplotnum)
% Input:
%   mat: n channels by time
%   subplotnum: where in the subplot
%   figHandle: figure handle
% Output:
%   pow: n channels by frequency
%   

% delete rows with all NaN
dataMat = mat(all(~isnan(mat),2),:);

colMat = dataMat'; % all computation are done for each column
n  = size(colMat,1); % segment length
Fs = 1000; % sample rate
freq = Fs*(0:(n/2))/n;
% compute FFT for average cross correlation
hanwin = repmat(hanning(n),1,size(colMat,2)); % make hanning window (nChn x nT)
[foi, tickLoc, tickLabel,~,~] = getFoiLabel(2,128,150,2); % lowFreq, highFreq, numFreqs, linORlog)

% For matrices, the FFT operation is applied to each column.
tmpFFT = fft(detrend(colMat).*hanwin); % detrend removes the trend from each column of the matrix
tmpFFT = tmpFFT(1:1025,:); % get half of the frequencies
tmpPow = abs(tmpFFT)'.^2; % compute power, nChn x nFreq
for i = 1:size(tmpPow,1)
    pow(i,:) = interp1(freq,tmpPow(i,:),foi,'spline');
end
figure(figHandle)
subplot(numRow, numCol, subplotnum)
mn = nanmean(pow,1);
er = nanstd(pow,[],1)/sqrt(size(pow,1));
shadedErrorBar(1:numel(foi),mn,er);
%shadedErrorBar(1:numel(foi),pow2db(mn),pow2db(er)); % dB looks flattened

%hold on
%p1 = plot(intF,mn,'LineWidth',2);
% xlabel('Freq [Hz]')
% ylabel('Power [dB]')
set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)