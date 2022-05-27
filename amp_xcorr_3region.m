function [lags, crosscorr_12, max_crosscorr_lag_12,crosscorr_23, max_crosscorr_lag_23,crosscorr_13, max_crosscorr_lag_13]...
    =amp_xcorr_3region(eeg1,eeg2,eeg3,samp_freq,low_freq,high_freq,maxlag)
% amp_crosscorr filters two eeg signals between a specified frequency band,
% calculates the crosscorrelation of the amplitude envelope of the filtered signals
% and returns the crosscorrelation as an output.
% USAGE: [lags, crosscorr, max_crosscorr_lag]=amp_crosscorr(eeg1,eeg2,samp_freq,low_freq,high_freq)
%INPUTS:
% eeg1-vector containing local field potential from brain area 1
% eeg2-vector containing local field potential from brain area 2
% samp_freq-sampling frequency, in Hz, of eeg1 and eeg2
% low_freq-low cut off, in Hz, of the band pass filter that will be applied to eeg1 and eeg2
% high_freq-high cut off, in Hz, of the band pass filter that will be applied to eeg1 and eeg2
%OUTPUTS:
% lags-vector contaning lags from -100 ms to +100 ms, over which the
% crosscorrelation was done
% crosscorr-vector with the crosscorrelation of the amplitude of eeg1 eeg2
% after being filtered between low_freq and high_freq
% max_crosscorr_lag-lag at which the crosscorrelation peaks. 
% Negative max_crosscorr_lag indicates that eeg1 is leading eeg2.
% A phase of zero refers to the trough of the (eg.theta) cycle.

% 4/30/2018 downloaded from Josh Gordon's 2010 JNeuroscience Methods paper
% 5/2/2018 modified by Angel Huang: changed filter, added plot

% check inputs
if nargin ~=7
error('ERROR in amp_crosscorr. There must be 7 inputs. - USAGE: [lags, crosscorr, max_crosscorr_lag]= amp_crosscorr(eeg1,eeg2,eeg3,samp_freq,low_freq,high_freq,maxLag);')
end
% if nargout ~=3
% error('ERROR in amp_crosscorr. There must be 3 outputs. - USAGE: [lags, crosscorr, max_crosscorr_lag]=amp_crosscorr(eeg1,eeg2,eeg3,samp_freq,low_freq,high_freq,maxLag);')
% end
%check consistency of data
if length(eeg1)~= length(eeg2) || length(eeg2)~= length(eeg3)
error('ERROR in amp_crosscorr. eeg1 and eeg2 must be vectors of the same size;')
end
s=size(eeg1);
if min(s)~=1
error('ERROR in amp_crosscorr. signal must be one-dimensional vectors')
end
s=size(eeg2);
if min(s)~=1
error('ERROR in amp_crosscorr. signal must be one-dimensional vectors')
end
s=size(eeg3);
if min(s)~=1
error('ERROR in amp_crosscorr. signal must be one-dimensional vectors')
end

% higher order filter has sharper edges but requires signal to be > 3*order
% long; eg. order=1000 needs 3000 data points 
order = 2*round(samp_freq); %determines the order of the filter used, data length in sec/3
if mod(order,2)~= 0
order = order-1;
end


Nyquist=floor(samp_freq/2);%determines nyquist frequency

% design filter (same result as fir1 function)
dbp = designfilt('bandpassfir', 'FilterOrder', order, ...
            'CutoffFrequency1',low_freq,'CutoffFrequency2',high_freq, ...
            'SampleRate', samp_freq);

filtered1 = filtfilt(dbp,eeg1); %filters eeg1 between low_freq and high_freq
filtered2 = filtfilt(dbp,eeg2); %filters eeg2 between low_freq and high_freq
filtered3 = filtfilt(dbp,eeg3); %filters eeg2 between low_freq and high_freq

filt_hilb1 = hilbert(filtered1); %calculates the Hilbert transform of eeg1
amp1_dc = abs(filt_hilb1);%calculates the instantaneous amplitude of eeg1 filtered between low_freq and high_freq
amp1=amp1_dc-mean(amp1_dc); %removes mean of the signal because the DC component of a signal does not change the correlation
filt_hilb2 = hilbert(filtered2);%calculates the Hilbert transform of eeg2
amp2_dc = abs(filt_hilb2);%calculates the instantaneous amplitude of eeg2 filtered between low_freq and high_freq
amp2=amp2_dc-mean(amp2_dc);
filt_hilb3 = hilbert(filtered3);%calculates the Hilbert transform of eeg2
amp3_dc = abs(filt_hilb3);%calculates the instantaneous amplitude of eeg2 filtered between low_freq and high_freq
amp3=amp3_dc-mean(amp3_dc);

[crosscorr_12,lags]=xcorr(amp1, amp2,maxlag,'coeff'); %calculates crosscorrelations between amplitude vectors
% amp1 and amp2 must have the same length to use 'coeff'
% maxlag = round(samp_freq/10): limits the lag range to be [–maxlag, maxlag]. 
% eg. if samp_freq = 1000, this will calculate lag from -100 to 100ms.
% 'coeff': Normalizes the sequence so that the autocorrelations at zero lag equal 1
lags=(lags./samp_freq)*1000; %converts lags to miliseconds, [-100,-99...99,100], always the same
g_12=find(crosscorr_12==max(crosscorr_12));%identifies index where the crosscorrelation peaks
max_crosscorr_lag_12=lags(g_12);%identifies the lag at which the crosscorrelation peaks

[crosscorr_23,lags]=xcorr(amp2, amp3,maxlag,'coeff'); %calculates crosscorrelations between amplitude vectors
lags=(lags./samp_freq)*1000; %converts lags to miliseconds
g_23=find(crosscorr_23==max(crosscorr_23));%identifies index where the crosscorrelation peaks
max_crosscorr_lag_23=lags(g_23);%identifies the lag at which the crosscorrelation peaks

[crosscorr_13,lags]=xcorr(amp1, amp3,maxlag,'coeff'); %calculates crosscorrelations between amplitude vectors
lags=(lags./samp_freq)*1000; %converts lags to miliseconds
g_13=find(crosscorr_13==max(crosscorr_13));%identifies index where the crosscorrelation peaks
max_crosscorr_lag_13=lags(g_13);%identifies the lag at which the crosscorrelation peaks


% plot raw LFP, filtered LFP and envelope
doPlot = 0;
if doPlot == 1
    screensize = get( groot, 'Screensize' );
    fig = figure('Position',[10 50 screensize(3)-100 screensize(4)-150]);
    
    red_rgb = [1 0 0]; %PPC
    blue_rgb = [0 0 1]; %LPl
    green_rgb = [0 1 0]; %VC
    
    tvec = 0:1/samp_freq:(length(eeg1)-1)/samp_freq;
    
    subplot(2, 3, 1); 
        plot(tvec, eeg1, 'color', 0.6*[1 1 1]); title('PPC'); hold on; 
        plot(tvec, filtered1,'linewidth', 2, 'color', red_rgb); 
        plot(tvec, amp1_dc, 'linewidth', 2, 'color', 0.6*red_rgb); 
        legend('raw LFP', 'filt', 'envelope'); xlabel('time [s]'); ylabel('uV'); xlim([tvec(1) tvec(end)]);
    subplot(2, 3, 2); 
        plot(tvec, eeg2,  'color', 0.6*[1 1 1]); title('LPl'); hold on; 
        plot(tvec, filtered2,'linewidth', 2, 'color', blue_rgb);  
        plot(tvec, amp2_dc, 'linewidth', 2, 'color', 0.6*blue_rgb); 
        legend('raw LFP', 'filt', 'envelope'); xlabel('time [s]'); ylabel('uV'); xlim([tvec(1) tvec(end)]);
    subplot(2, 3, 3); 
        plot(tvec, eeg3,  'color', 0.6*[1 1 1]); title('VC'); hold on; 
        plot(tvec, filtered3,'linewidth', 2, 'color', green_rgb);  
        plot(tvec, amp3_dc, 'linewidth', 2, 'color', 0.6*green_rgb); 
        legend('raw LFP', 'filt', 'envelope'); xlabel('time [s]'); ylabel('uV'); xlim([tvec(1) tvec(end)]);
    subplot(2, 3, 4); 
        plot(tvec, amp1_dc,  'color', red_rgb,'linewidth', 2); title('Amplitude'); hold on; 
        plot(tvec, amp2_dc,  'color', blue_rgb,'linewidth', 2); 
        plot(tvec, amp3_dc,  'color', green_rgb,'linewidth', 2); 
        legend('PPC','LPl','V1'); ylabel('uV'); xlim([tvec(1) tvec(end)]);
    subplot(2,3,5);        
        plot(tvec, amp1,  'color', red_rgb,'linewidth', 2); title('Amplitude'); hold on; 
        plot(tvec, amp2,  'color', blue_rgb,'linewidth', 2); 
        plot(tvec, amp3,  'color', green_rgb,'linewidth', 2); 
        legend('PPC','LPl','V1'); ylabel('uV'); xlim([tvec(1) tvec(end)]);    
    
    subplot(2, 3, 6);
        plot(lags, crosscorr_12,'color','r','linewidth',2),hold on %plots crosscorrelations
        plot(lags, crosscorr_23,'color','b','linewidth',2),hold on %plots crosscorrelations
        plot(lags, crosscorr_13,'color','g','linewidth',2),hold on %plots crosscorrelations
        vline(max_crosscorr_lag_12, 'r:');
        vline(max_crosscorr_lag_23, 'b:');
        vline(max_crosscorr_lag_13, 'g:');
        vline(0, 'k--');
        plot(lags(g_12),crosscorr_12(g_12),'r^','markerfacecolor','r','markersize',10)%plots marker at the peak of the cross correlation
        plot(lags(g_23),crosscorr_23(g_23),'b^','markerfacecolor','b','markersize',10)%plots marker at the peak of the cross correlation
        plot(lags(g_13),crosscorr_13(g_13),'g^','markerfacecolor','g','markersize',10)%plots marker at the peak of the cross correlation

        %plot([0 0],[1.05*max(crosscorr) 0.95*min(crosscorr)],'color',[0 0 0],'linestyle',':', 'linewidth',2) %plots dashed line at zero lag
        %set(gca,'xtick',[-100 -50 0 50 100])
        axis tight; box off; %xlim([-101 100])
        xlabel('Lag [ms]','fontsize',14)
        ylabel('Crosscorrelation','fontsize',14)
        legend('PPC-LPl lag','LPl-V1 lag','PPC-V1 lag');
        title('Theta amplitude crosscorrelation');
        
  savefig(fig, 'amp_xcorr_3region_7sec_0.fig','compact');
  saveas(fig, 'amp_xcorr_3region_7sec_0.png');      
end
