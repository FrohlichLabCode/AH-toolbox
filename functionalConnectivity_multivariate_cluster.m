function funcCon = EN_functionalConnectivity_multivariate_cluster(xser,yser, zser,fs,event,twin,chn1,chn2, chn3, saveRootPath)
% This function computes a various metrics of functional and effective
% connectivity time-locked to certain events, including:
% 	* Spectral power
% 	* Phase locking value
% 	* Coherence
% 	* Imaginary coherence
% 	* Phase slope index
% 	* Granger causality
% 
% This code also generates a single plot where you can compare all of the
% above-mentioned metrics. 
% Input: - xser  (time series from one brain region)
%        - yser  (time series from another brain region)
%        - fs    (sample rate of time series)
%        - event (vector of event times in seconds)
%        - twin  (time window for spectral analysis, eg twin = [-2 5])
% Output - funcCon (structure containing all information on functional/effective connectivity analyses)
%
% For this script you will need to have the MVGC toolbox in your path. You
% can find the website for the toolbox here: http://users.sussex.ac.uk/~lionelb/MVGC/
% Or you can also just copy the code from the following folder in my
% codebase: C:\Users\FrohlichLab\Dropbox (Frohlich Lab)\Codebase\CodeIain\mainCode\Code\mvgc_v1.0
% I.S. 2017
% E.N. 10/27/2017 Modified for multivariate Granger Causality
% A.H. 4/25/2018 Added time-inverted multivariate Granger Causality

% initialize data structure
funcCon = struct;

% Define frequencies of interest. Linear spacing for Phase slope index, and
% logarithmic spacing for all other methods.
numFreqs = 100;
lowFreq  = 2;
highFreq = 30;
foi      = linspace(lowFreq,highFreq,numFreqs); % linear spacing
% foi   = logspace(log10(lowFreq),log10(highFreq),numFreqs); % log spacing

% reject noise in recording (high amplitude noise can destroy coherence estimates in particular)
xsig = sub_rejectNoise(xser,fs,2,1); 
ysig = sub_rejectNoise(yser,fs,2,1); 
zsig = sub_rejectNoise(zser,fs,2,1); 

% Subtract event related potential to look at only induced oscillations
subtractERP = 1;
if subtractERP == 1
    % Finish this part next week... 
end


% % % % % % % % % % % % % % % % % % % % % % % % % %
% Compute pairwise-conditional Granger Causality  %
% % % % % % % % % % % % % % % % % % % % % % % % % %

funcCon = sub_grangerCausality(funcCon,xser,yser,zser,isnan(xsig),isnan(ysig),...
    isnan(zsig),event,fs,twin,foi, saveRootPath);

tvecGC = funcCon.grangerCausality.tvec;
doPlot = 0;
if doPlot == 1
    %     screensize = get( groot, 'Screensize' );
    fig = figure('Position',[-1908          38        1438         930]);
    % Compute ticks for plotting
    fois = [0.5 1 2 4 8 16 32 64 128];
    for fi = 1:numel(fois)
        [bi,bb] = sort(abs(foi-fois(fi)));
        tickLoc(fi) = bb(1);
    end
    tickLabel = {'0.5','1','2','4','8','16','32','64','128'};
    
    
    
    
    myclim = 0.4;
    myylim = 75;
    subplot(231);   imagesc(tvecGC,foi,real(funcCon.grangerCausality.X_to_Z));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: X to Y')
    set(gca,'YDir','normal','TickDir','out')
    cl = colorbar('northoutside'); ylabel(cl,'GC: PPC to VC','FontSize',12)
    caxis([0 myclim]);    ylim([30 myylim]); % 90+ Hz has saturated values
    
    subplot(232);   imagesc(tvecGC,foi,real(funcCon.grangerCausality.Z_to_X));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: X to Y')
    set(gca,'YDir','normal','TickDir','out')
    cl = colorbar('northoutside'); ylabel(cl,'GC:  VC to PPC','FontSize',12)
    caxis([0 myclim]);    ylim([30 myylim]); % 90+ Hz has saturated values
    
    subplot(233);   imagesc(tvecGC,foi,real(funcCon.grangerCausality.X_to_Y));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: X to Y')
    set(gca,'YDir','normal','TickDir','out')
    cl = colorbar('northoutside'); ylabel(cl,'GC: PPC to LP','FontSize',12)
    caxis([0 myclim]);    ylim([30 myylim]); % 90+ Hz has saturated values
    
    subplot(234);   imagesc(tvecGC,foi,real(funcCon.grangerCausality.Y_to_X));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: X to Y')
    set(gca,'YDir','normal','TickDir','out')
    cl = colorbar('northoutside'); ylabel(cl,'GC: LP to PPC','FontSize',12)
    caxis([0 myclim]);    ylim([30 myylim]); % 90+ Hz has saturated values
    
    subplot(235);   imagesc(tvecGC,foi,real(funcCon.grangerCausality.Y_to_Z));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: X to Y')
    set(gca,'YDir','normal','TickDir','out')
    cl = colorbar('northoutside'); ylabel(cl,'GC: LP to VC','FontSize',12)
    caxis([0 myclim]);    ylim([30 myylim]); % 90+ Hz has saturated values
    
    subplot(236);   imagesc(tvecGC,foi,real(funcCon.grangerCausality.Z_to_Y));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: X to Y')
    set(gca,'YDir','normal','TickDir','out')
    cl = colorbar('northoutside'); ylabel(cl,'GC: VC to LP','FontSize',12)
    caxis([0 myclim]);    ylim([30 myylim]); % 90+ Hz has saturated values
    colormap(jet)
    
        if ~exist([saveRootPath 'figures'],'dir'); mkdir([saveRootPath 'figures']); end
        saveas(fig,[saveRootPath, 'figures/' num2str(chn1) '-' num2str(chn2) '-' num2str(chn3) '.png']);
    
end

return

function wav = sub_makeWavelet(foi,Fs)
% This function generates complex morlet wavelets that are Gaussian shaped
% in both the time and frequency domains. The equtions used to generate
% wavelets are taken from Tallon-Baundry et al (1996) J Neuroscience.
% 
% Inputs:  foi - vector of center frequencies to generate wavelets
%          Fs  - sample frequency of signal that is to be convolved with wavelets
% Outputs: wav - a cell array that contains the complex morlet wavelets 
% I.S. 2016 

q  = 7; % Wavelet width in cycles
w  = 3; % Width of wavelet taper in terms of standard deviation

wav = cell(numel(foi),1);
for f = 1:numel(foi)
    sf     = foi(f)/q;                   % standard deviation in the frequency domain
    st     = 1/(2*pi*sf);                % standard deviation in the time domain
    t      = -w*st:1/Fs:w*st;            % time vector for wavelet calculation
    A      = 1/sqrt(st*sqrt(pi));        % normalisation factor
    tap    = (A*exp(-t.^2/(2*st^2)));    % Gaussian window
    wav{f} = tap.*exp(2*1i*pi*foi(f)*t); % wavelet formula
    E      = sum(abs(wav{f}).^2);        % energy of wavelet
    wav{f} = wav{f}./sqrt(E);            % normalise wavelet energy to 1
end
return

function sigOut = sub_rejectNoise(sigIn,fs,winLen,rejThresh)
sigOut      = sigIn;
% delWin      = ones(1,round(fs*winLen));                   % window to cut out
% delInd      = (abs(zscore(sigIn)) > rejThresh);           % Samples that surpass threshold
% delVec      = (conv(double(delInd),delWin,'same') > 0);   % Convolve to smooth window and detect non-zero samples
% sigOut(delVec) = nan;                                     % Noisy samples change to NaN's

rejThresh = 200; % 500uV threshold for rejection
[b,a] = butter(4,40/(fs/2),'high'); % define highpass filter at 40Hz
hfSig = filtfilt(b,a,sigIn); % filtyer signal
delWin      = ones(1,round(fs*winLen));                   % window to cut out
delInd      = (abs((hfSig)) > rejThresh);           % Samples that surpass threshold
delVec      = (conv(double(delInd),delWin,'same') > 0);   % Convolve to smooth window and detect non-zero samples
sigOut(delVec) = nan;

return

function funcCon = sub_grangerCausality(funcCon,xsig,ysig,zsig,xnan,ynan,znan,event,fs,twin,foi, saveRootPath)

% set up priors
regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)
morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 20;     % maximum model order for model order estimation
acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)
fres      = [];     % frequency resolution (empty for automatic calculation)
segLength = 1;      % window length for computation of Granger Causality
newFs     = 200; % This is very important for GC! 

% define low pass filter at 100Hz
nyqFreq = fs/2;
[b,a]   = butter(2,100/nyqFreq,'low');
xfilt   = filtfilt(b,a,xsig);
yfilt   = filtfilt(b,a,ysig);
zfilt   = filtfilt(b,a,zsig);

% resample data to sample rate of newFs (defined above)
idat = resample(xfilt,newFs,fs);
jdat = resample(yfilt,newFs,fs);
zdat = resample(zfilt,newFs,fs);
inan = round(resample(double(xnan),newFs,fs));
jnan = round(resample(double(ynan),newFs,fs));
znan = round(resample(double(znan),newFs,fs));

% to get time resolution of GC, divide window of interests into 1 sec
% sliding windows, and calculate conditional GC for each sliding window 
stepSize  = 0.1; % sliding window increments
stepCen   = twin(1):stepSize:twin(2); % sliding window centers for each GC calculation
recLength = numel(idat)/newFs; % whole recording length in sec
halfWinSamp = (segLength*newFs)/2; 

X2Y = nan(numel(foi),numel(stepCen));
Y2X = nan(numel(foi),numel(stepCen));
X2Z = nan(numel(foi),numel(stepCen));
Z2X = nan(numel(foi),numel(stepCen));
Y2Z = nan(numel(foi),numel(stepCen));
Z2Y = nan(numel(foi),numel(stepCen));

X2Y_invert = nan(numel(foi),numel(stepCen));
Y2X_invert = nan(numel(foi),numel(stepCen));
X2Z_invert = nan(numel(foi),numel(stepCen));
Z2X_invert = nan(numel(foi),numel(stepCen));
Y2Z_invert = nan(numel(foi),numel(stepCen));
Z2Y_invert = nan(numel(foi),numel(stepCen));


% imat = 0; % for parfor
% jmat = 0; % for parfor
% zmat = 0; % for parfor
% X = 0;
for istep = 1:numel(stepCen) % loop through each sliding time window centers to calculate GC
    c = 0; % counting variable
    clear imat jmat zmat % due to parfor
    % fill up matrices with data 
    for iev = 1:numel(event) % loop through each trial
        tsamps = round(twin*fs);
        % skip if we have window range issues
        if event(iev) < abs(twin(1))+segLength; continue; end
        if event(iev) > recLength-twin(2)-segLength; continue; end
        samp      = round((stepCen(istep)+event(iev))*newFs);
        
        itmp = inan(samp-halfWinSamp:samp+halfWinSamp);
        jtmp = jnan(samp-halfWinSamp:samp+halfWinSamp);
        ztmp = znan(samp-halfWinSamp:samp+halfWinSamp);
        if sum(itmp) + sum(jtmp) + sum(ztmp) == 0 % only use data that have no noise (ie, no nan's)
            c = c + 1; % count for trial
            imat(:,c) = idat(samp-halfWinSamp:samp+halfWinSamp);
            jmat(:,c) = jdat(samp-halfWinSamp:samp+halfWinSamp);
            zmat(:,c) = zdat(samp-halfWinSamp:samp+halfWinSamp);
        else
            continue
        end
    end
    clear X %due to parfor
    X(1,:,:)  = imat;  % region x time x trial
    X(2,:,:)  = jmat;
    X(3,:,:)  = zmat;
    numSeg    = c;
%     imat=0; jmat=0; zmat=0;
    
    %% Model order estimation (<mvgc_schema.html#3 |A2|>)

    % compute information criterion
    nvars = size(X,1);
    [AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode, 0);
    amo = 10; % actual model order
    
    % Select model order.
    if strcmpi(morder,'actual')
        morder = amo;
%         fprintf('\nusing actual model order = %d\n',morder);
    elseif strcmpi(morder,'AIC')
        morder = moAIC;
%         fprintf('\nusing AIC best model order = %d\n',morder);
    elseif strcmpi(morder,'BIC')
        morder = moBIC;
%         fprintf('\nusing BIC best model order = %d\n',morder);
    else
%         fprintf('\nusing specified model order = %d\n',morder);
    end
    
    %% VAR model estimation (<mvgc_schema.html#3 |A2|>)
    
    % Estimate VAR model of selected order from data.
    % Fit autoregressive model to data
    [A,SIG] = tsdata_to_var(X,morder,regmode);
%     X = 0;
    assert(~isbad(A),'VAR estimation failed'); % Check for failed regression
    % NOTE: at this point we have a model and are finished with the data! - all
    % subsequent calculations work from the estimated VAR parameters A and SIG.
    %% Autocovariance calculation (<mvgc_schema.html#3 |A5|>)
    
    % The autocovariance sequence drives many Granger causality calculations (see
    % next section). Now we calculate the autocovariance sequence G according to the
    % VAR model, to as many lags as it takes to decay to below the numerical
    % tolerance level, or to acmaxlags lags if specified (i.e. non-empty).
    % return autocovariance sequence
    [G,info] = var_to_autocov(A,SIG,acmaxlags);
    
    % var_info(info,true); % report results (and bail out on error) --skip, don't want to abort
    
    %% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)
    
    % Calculate time-domain pairwise-conditional causalities - this just requires
    % the autocovariance sequence.
    % compute Granger Causality based on autocovariance
    %f: a 3-dim numerical matrix with
    % first index representing target ("to"), second index source ("from")
    % quantities and third index frequencies - typically spectral causalities
    f = autocov_to_spwcgc(G,fres);
    try
        assert(~isbad(f,false),'spectral GC calculation failed');
        freqRes = size(f,3)-1;
        freqs   = linspace(0,newFs/2,freqRes+1)'; % compute frequencies
    
        % interpolate to frequencies of interest
        X2Y(:,istep) = interp1(freqs,squeeze(f(2,1,:)),foi,'spline');
        Y2X(:,istep) = interp1(freqs,squeeze(f(1,2,:)),foi,'spline');    
        X2Z(:,istep) = interp1(freqs,squeeze(f(3,1,:)),foi,'spline');    
        Z2X(:,istep) = interp1(freqs,squeeze(f(1,3,:)),foi,'spline');    
        Y2Z(:,istep) = interp1(freqs,squeeze(f(3,2,:)),foi,'spline');    
        Z2Y(:,istep) = interp1(freqs,squeeze(f(2,3,:)),foi,'spline');    
    catch
    end    

    % calculate invert time series Granger causality
    [A2,SIG2] = tsdata_to_var(flip(X,2),morder,regmode);
    assert(~isbad(A2),'VAR estimation failed'); % Check for failed regression
    [G2,info] = var_to_autocov(A2,SIG2,acmaxlags);
    f2 = autocov_to_spwcgc(G2,fres);
    try
        assert(~isbad(f2,false),'spectral GC calculation failed');
        freqRes = size(f2,3)-1;
        freqs   = linspace(0,newFs/2,freqRes+1)'; % compute frequencies
    
        % interpolate to frequencies of interest
        X2Y_invert(:,istep) = interp1(freqs,squeeze(f2(2,1,:)),foi,'spline');
        Y2X_invert(:,istep) = interp1(freqs,squeeze(f2(1,2,:)),foi,'spline');    
        X2Z_invert(:,istep) = interp1(freqs,squeeze(f2(3,1,:)),foi,'spline');    
        Z2X_invert(:,istep) = interp1(freqs,squeeze(f2(1,3,:)),foi,'spline');    
        Y2Z_invert(:,istep) = interp1(freqs,squeeze(f2(3,2,:)),foi,'spline');    
        Z2Y_invert(:,istep) = interp1(freqs,squeeze(f2(2,3,:)),foi,'spline');    
    catch
    end  
    
    
end

% try
    funcCon.grangerCausality.X_to_Y = X2Y;
    funcCon.grangerCausality.Y_to_X = Y2X;
    funcCon.grangerCausality.X_to_Z = X2Z;
    funcCon.grangerCausality.Z_to_X = Z2X;
    funcCon.grangerCausality.Y_to_Z = Y2Z;
    funcCon.grangerCausality.Z_to_Y = Z2Y;
    funcCon.grangerCausality.X_to_Y_invert = X2Y_invert;
    funcCon.grangerCausality.Y_to_X_invert = Y2X_invert;
    funcCon.grangerCausality.X_to_Z_invert = X2Z_invert;
    funcCon.grangerCausality.Z_to_X_invert = Z2X_invert;
    funcCon.grangerCausality.Y_to_Z_invert = Y2Z_invert;
    funcCon.grangerCausality.Z_to_Y_invert = Z2Y_invert;
    funcCon.grangerCausality.tvec   = stepCen;
%     GCtvec = stepCen;
%     save([saveRootPath, 'GCtvecxx.mat'], 'stepCen'); % I needed to do this to handle some errors about tvec with parfor
    fprintf('Functional connectivity of this pair is done!!!!!!!!!!!!!')
% catch
% end
return

testing = 0;
if testing == 1
    % Let's generate some surrogate data for testing this code
    numEvs = 100; % define number of events
    fs     = 1e3; % sample rate
    fband  = [30 60]; % frequency band of surrogate interaction
    evDur  = 2; % surrogate stimulus duration
    evs    = (1:numEvs)*10+rand(1,numEvs); % Define event times
    reclength = evs(end)+evDur+5; % length of vector to make (in seconds)
    recSamps  = round(reclength*fs); % number of samples in surrogate data vector
    [c,d]     = butter(2,0.5/(fs/2),'high'); % highpass filter to remove low freq components
    x         = filter(c,d,pinknoise(recSamps)); % highpass filter pink noise (1/f)
    y         = filter(c,d,pinknoise(recSamps)); % highpass filter pink noise (1/f)
    [b,a]     = butter(2,fband/(fs/2),'bandpass'); % bandpass filter for adding band limited signals
    s         = filter(b,a,randn(1,recSamps))*2; % surrogate band limited signal
    timeLag   = -round(0.005*fs); % time lag between two surrogate signals (in seconds)
    randSamps = round((rand(1,numEvs)*(reclength-10))*fs); % samples to take surrogate oscillatory signal

    % Loop through events and add the band-limited surrogate data
    for iev = 1:numEvs
        samp  = round(evs(iev)*fs);
        rsamp = randSamps(iev);
        x(samp:samp+(evDur*fs)) = x(samp:samp+(evDur*fs)) + s(rsamp:rsamp+(evDur*fs)); % add band limited data
        shiftSamp = rsamp + timeLag;
        y(samp:samp+(evDur*fs)) = y(samp:samp+(evDur*fs)) + s(shiftSamp:shiftSamp+(evDur*fs)) + rand(1,evDur*fs+1); % add band limited data with offset plus some noise
    end
end
%funcCon = is_functionalConnectivity(x,y,fs,evs,[-2 4])

