%mat = csvread('F:\DB\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeEhsan\tmpforFlavio\1.csv',2,1);
%mat = mat(:,1);
mat = xsig;
fs = 100;
%fs  = 1e6;%1MHz
tvec = -0.6:1/fs:(length(mat)-1)/fs-0.6;
window   = 0.9*fs; % frequency resolution, higher is higher resolution
noverlap = round(0.1*window);

%% define foi to compute
linORlog = 1;
if linORlog == 1
    numFreqs = 100;
    lowFreq  = 2;
    highFreq = 30;
    foi      = linspace(lowFreq,highFreq,numFreqs); % linear spacing
    fois = [2, 5:5:highFreq];
    tickLabel = string(fois); % generate a string array matches fois {"5","10"...}
elseif linORlog == 2
    numFreqs = 80;
    lowFreq  = 2;
    highFreq = 128;
    foi   = logspace(log10(lowFreq),log10(highFreq),numFreqs); % log spacing
    fois = 2.^(log2(lowFreq):1:log2(highFreq)); %[2 4 8 12 16 32 64 128];
    tickLabel = string(fois);
end
% [tickLoc, tickLabel] = getTickLabel(lowFreq, highFreq, numFreqs, linORlog);
sinW = max(mat)*sin(2*pi*10*tvec);
%% plot
fig = figure('name','Example signal');
subplot(4,1,1)
  plot(tvec,mat); hold on
  plot(tvec, sinW, 'color', 'r', 'linestyle', '--')
  xlabel('Time [sec]');
  title('Raw trace');
  legend('signal', 'sine wave')
subplot(4,1,2) 
%   nfft=100; 
  %window = hanning(nfft); %noverlap = nfft/2;
  [s,f,t] = spectrogram(xsig,window,noverlap,fs);
  
  
  [sW,f,t] = spectrogram(sinW,window,noverlap,[],fs);
  power = abs(s).^2;
  powerW = abs(sW).^2;
  imagesc(tvec,foi,pow2db(power));
  cl = colorbar; ylabel(cl,['Power [dB]']);
  ylabel('Frequency [Hz]');
  xlabel('Time [sec]');
  title('Power spectrogram');
  xlim([-0.6,0.6]);
  caxis([80,120]);
  set(gca,'YDir','normal');%,'TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel);
subplot(4,1,3)
  plot(foi,nanmean(power,2), 'linewidth', 2); hold on
  plot(foi,nanmean(powerW,2), 'color', 'r', 'linestyle', '--', 'linewidth', 2);
  ylabel('Power');
  xlabel('Frequency [Hz]');
  title('Power spectrum');
  xlim([2 15 ])
subplot(4,1,4)
  plot(foi,10*log10(nanmean(power,2)), 'linewidth', 2); hold on
  plot(foi,10*log10(nanmean(powerW,2)), 'color', 'r', 'linestyle', '--', 'linewidth', 2);
  ylabel('Power [dB]');
  xlabel('Frequency [Hz]');
  title('Power spectrum');
  xlim([2 15 ])
fig.Color = 'white';