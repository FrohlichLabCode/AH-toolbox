x = readtable('C:\Users\Angel\Desktop\signal_1_downsampled.dat');
y = readtable('C:\Users\Angel\Desktop\signal_2_downsampled.dat');
xsig = table2array(x);
ysig = table2array(y);
fs = 100;
linORlog = 2;
lowFreq = 2;
highFreq = 50;
numFreqs = 50;
% Define frequencies of interest. Linear spacing for Phase slope index, and spacing for all other methods.
if     linORlog == 1
        foi      = linspace(lowFreq,highFreq,numFreqs); % linear spacing
elseif linORlog == 2
        foi      = logspace(log10(lowFreq),log10(highFreq),numFreqs); % log spacing
end
morWav = is_makeWavelet(foi,fs);
for f = 1:numFreqs
    xspec(f,:) = conv(xsig,morWav{f},'same');
    yspec(f,:) = conv(ysig,morWav{f},'same');
end
xmat = xspec;
ymat = yspec;

[s,f,t] = spectrogram(mat,window,noverlap,fs);

Sxy = xmat.*conj(ymat); % multiply x by complex conjugate of y
Sxx = xmat.*conj(xmat); % multiply x by complex conjugate of x
Syy = ymat.*conj(ymat); % multiply y by complex conjugate of y
Cy  = Sxy./(sqrt(Sxx.*Syy)); % coherency formula
funcCon.coherency          = Cy; % might want to keep the complex part to look at phase lags 
funcCon.coherence          = abs(Cy);

figure()
plot(1:numel(foi),funcCon.coherence(:,1));
xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Coherence')
ylim([tickLoc(1) tickLoc(end)]);
%caxis([0 0.8]);
set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel);
cl = colorbar('northoutside'); ylabel(cl,'Coherence','FontSize',12);

% plv
phaseDiff   = angle(xmat) - angle(ymat); % compute phase angle lag between signals (note: not the diff between analytical signals)
plv     = abs(nanmean(exp(1i*phaseDiff),2)); % PLV formula
plot(plv);
figure()
subplot(411);plot(angle(xmat(:,25)));title('xPhase')
subplot(412);plot(angle(ymat(:,25)));title('yPhase')
subplot(413);plot(phaseDiff(:,25));title('PhaseDiff')
subplot(414);plot(plv);title('plv')
