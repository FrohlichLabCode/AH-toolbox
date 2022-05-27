

animals = {'0114','0116','0124','0125'};
pathDir       = 'E:\';
lowFreq      = 0.5;
highFreq     = 128;
numFreqs     = 80;
foi          = logspace(log10(lowFreq),log10(highFreq),numFreqs);
logF         = logspace(log10(lowFreq),log10(highFreq),1000);
newFs        = 200;

% set up priors
regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)
morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 20;     % maximum model order for model order estimation
acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)
tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')
seed      = 0;      % random seed (0 for unseeded)
fres      = [];     % frequency resolution (empty for automatic calculation)
segLength = 1;
numPerm   = 100;    % number of random permutations used to calculate significance
bsize     = [];
nperm     = 1000;

count = 0;
for ianimal = 1:numel(animals)
    animalCode = animals{ianimal};
    socialDarkDir = dir([pathDir animalCode '\' animalCode '_imagesVideo*']);
    
    for irec = 1:numel(socialDarkDir);
        
        recName = socialDarkDir(irec).name;
        recPath = [pathDir animalCode '\' recName '\'];
        load([recPath 'triggerData'])
        
        try
            % read in the log file
            fileName = ['E:\Presentation_LogFiles\' recName(1:end-7) '.log'];
            fileID   = fopen(fileName);
            formatSpec = '%f %s %s %f %f %f %f %f %f %f %s %f';
            LogFile = textscan(fileID,formatSpec,'HeaderLines',5,'Delimiter', '\t');
            fclose(fileID);
        catch; continue; end
        
        display(recName)
        count = count + 1;
        load([recPath 'lfp\lfpMat'])
        keepChans = keepPulChans(animalCode);
        if size(lfpMat,1) >48; pulChans = keepChans + 40; else pulChans = keepChans + 32; end
        ppcChans = 1:32;
        
        % detect images
        cellPic = cellfun(@regexp,LogFile(2),{'Picture'},'UniformOutput',false);
        picInd  = find(cellfun(@isempty,cellPic{:}) == 0);
        % detect videos
        cellVid = cellfun(@regexp,LogFile(2),{'Video'},'UniformOutput',false);
        vidInd  = find(cellfun(@isempty,cellVid{:}) == 0);
        
        % detect stim onsets
        onset  = find(diff(triggerData(1,:))==1);
        offset = find(diff(triggerData(2,:))==1);
        
        picOnset = (onset(picInd)/Fs);
        vidOnset = (onset(vidInd)/Fs);
        
        [b,a] = butter(2,100/(lfpFs/2),'low');
        fdat  = filtfilt(b,a,lfpMat')';
        
        idat = resample(median(fdat(ppcChans,:)),newFs,lfpFs);
        jdat = resample(median(fdat(pulChans,:)),newFs,lfpFs);
        
        win       = [-5 15];
        segLength = 1;
        stepSize  = 0.25;
        stepCen   = win(1):stepSize:win(2);
        recLength = numel(idat)/newFs;
        halfWinSamp = (newFs*segLength)/2;
        nobs      = round(segLength*newFs);
        
        PPC2Pul = nan(numel(stepCen),numel(logF));
        Pul2PPC = nan(numel(stepCen),numel(logF));
        for istep = 1:numel(stepCen)
            c         = 0;
            clear imat jmat
            for isac = 1:numel(picOnset)
                if picOnset(isac) < 5.5; continue; end
                if picOnset(isac) > recLength-15.5; continue; end
                c = c + 1;
                samp      = round((picOnset(isac)+stepCen(istep))*newFs);
                imat(:,c) = idat(samp-halfWinSamp:samp+halfWinSamp);
                jmat(:,c) = jdat(samp-halfWinSamp:samp+halfWinSamp);
            end
            clear X
            X(1,:,:)  = imat;
            X(2,:,:)  = jmat;
            numSeg    = c;
            
            nvars = size(X,1);
            [AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
            amo = 10; % actual model order
            
            % Select model order.
            if     strcmpi(morder,'actual')
                morder = amo;
                fprintf('\nusing actual model order = %d\n',morder);
            elseif strcmpi(morder,'AIC')
                morder = moAIC;
                fprintf('\nusing AIC best model order = %d\n',morder);
            elseif strcmpi(morder,'BIC')
                morder = moBIC;
                fprintf('\nusing BIC best model order = %d\n',morder);
            else
                fprintf('\nusing specified model order = %d\n',morder);
            end
            
            % morder = 10;
            
            [A,SIG] = tsdata_to_var(X,morder,regmode);
            assert(~isbad(A),'VAR estimation failed');
            
            [G,info] = var_to_autocov(A,SIG,acmaxlags);
            
            try var_info(info,true); % report results (and bail out on error)
            catch; continue; end
            
            F = autocov_to_pwcgc(G);
            assert(~isbad(F,false),'GC calculation failed');
            
            pval = mvgc_pval(F,morder,nobs,numSeg,1,1,nvars-2,tstat); % take careful note of arguments!
            sig  = significance(pval,alpha,mhtc);
            
            f = autocov_to_spwcgc(G,fres);
            assert(~isbad(f,false),'spectral GC calculation failed');
            freqRes        = size(f,3)-1;
            freqs          = linspace(0,newFs/2,freqRes+1)';
            
            PPC2Pul(istep,:) = interp1(freqs,squeeze(f(2,1,:)),logF,'spline');
            Pul2PPC(istep,:) = interp1(freqs,squeeze(f(1,2,:)),logF,'spline');
        end
        pop_PPC2Pul_pic(count,:,:) = PPC2Pul;
        pop_Pul2PPC_pic(count,:,:) = Pul2PPC;
        
        
        PPC2Pul = nan(numel(stepCen),numel(logF));
        Pul2PPC = nan(numel(stepCen),numel(logF));
        for istep = 1:numel(stepCen)
            c         = 0;
            clear imat jmat
            for isac = 1:numel(vidOnset)
                if vidOnset(isac) < 5.5; continue; end
                if vidOnset(isac) > recLength-15.5; continue; end
                c = c + 1;
                samp      = round((vidOnset(isac)+stepCen(istep))*newFs);
                imat(:,c) = idat(samp-halfWinSamp:samp+halfWinSamp);
                jmat(:,c) = jdat(samp-halfWinSamp:samp+halfWinSamp);
            end
            clear X
            X(1,:,:)  = imat;
            X(2,:,:)  = jmat;
            numSeg    = c;
            
            nvars = size(X,1);
            [AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
            amo = 10; % actual model order
            
            % Select model order.
            if     strcmpi(morder,'actual')
                morder = amo;
                fprintf('\nusing actual model order = %d\n',morder);
            elseif strcmpi(morder,'AIC')
                morder = moAIC;
                fprintf('\nusing AIC best model order = %d\n',morder);
            elseif strcmpi(morder,'BIC')
                morder = moBIC;
                fprintf('\nusing BIC best model order = %d\n',morder);
            else
                fprintf('\nusing specified model order = %d\n',morder);
            end
            
            % morder = 10;
            
            [A,SIG] = tsdata_to_var(X,morder,regmode);
            assert(~isbad(A),'VAR estimation failed');
            
            [G,info] = var_to_autocov(A,SIG,acmaxlags);
            
            try var_info(info,true); % report results (and bail out on error)
            catch; continue; end
            
            F = autocov_to_pwcgc(G);
            assert(~isbad(F,false),'GC calculation failed');
            
            pval = mvgc_pval(F,morder,nobs,numSeg,1,1,nvars-2,tstat); % take careful note of arguments!
            sig  = significance(pval,alpha,mhtc);
            
            f = autocov_to_spwcgc(G,fres);
            assert(~isbad(f,false),'spectral GC calculation failed');
            freqRes        = size(f,3)-1;
            freqs          = linspace(0,newFs/2,freqRes+1)';
            
            PPC2Pul(istep,:) = interp1(freqs,squeeze(f(2,1,:)),logF,'spline');
            Pul2PPC(istep,:) = interp1(freqs,squeeze(f(1,2,:)),logF,'spline');
        end
        pop_PPC2Pul_vid(count,:,:) = PPC2Pul;
        pop_Pul2PPC_vid(count,:,:) = Pul2PPC;
    end
end

pop_PPC2Pul_pic(6,:,:) = nan;
pop_Pul2PPC_pic(6,:,:) = nan;
pop_PPC2Pul_vid(6,:,:) = nan;
pop_Pul2PPC_vid(6,:,:) = nan;

fois = [0.5 1 2 4 8 16 32 64 128];
for fi = 1:numel(fois)
    [bi,bb] = sort(abs(logF-fois(fi)));
    tickLoc(fi) = bb(1);
end

for irec = 1:size(pop_PPC2Pul_vid,1)
    for f = 1:size(pop_PPC2Pul_vid,3)
        vid_PPC2Pul(irec,:,f) = smooth(squeeze(pop_PPC2Pul_vid(irec,:,f)),3);
        vid_Pul2PPC(irec,:,f) = smooth(squeeze(pop_Pul2PPC_vid(irec,:,f)),3);
        pic_PPC2Pul(irec,:,f) = smooth(squeeze(pop_PPC2Pul_pic(irec,:,f)),3);
        pic_Pul2PPC(irec,:,f) = smooth(squeeze(pop_Pul2PPC_pic(irec,:,f)),3);
    end
end

mnPPC2Pul = flipud(rot90(squeeze(nanmean(pop_PPC2Pul_pic))));
mnPul2PPC = flipud(rot90(squeeze(nanmean(pop_Pul2PPC_pic))));

figure;
subplot(1,2,1)
% PPC 2 Pulvinar
imagesc(stepCen,1:numel(logF),mnPPC2Pul);
set(gca,'TickDir','out','YTick',tickLoc,'YDir','normal')
ax = gca;
ax.YTickLabel = {'0.5','1','2','4','8','16','32','64','128'};
xlabel('Time to saccade (s)')
ylabel('Frequency')
colormap(awesomemap2)
caxis([0 0.07])
ylim([200 800])
title('PPC to Pul image')

% Pulvinar 2 PPC
subplot(1,2,2)
% PPC 2 Pulvinar
imagesc(stepCen,1:numel(logF),mnPul2PPC);
set(gca,'TickDir','out','YTick',tickLoc,'YDir','normal')
ax = gca;
ax.YTickLabel = {'0.5','1','2','4','8','16','32','64','128'};
xlabel('Time to saccade (s)')
ylabel('Frequency')
colormap(awesomemap2)
caxis([0 0.1])
ylim([200 800])
title('Pul to PPC image')

mnPPC2Pul = flipud(rot90(squeeze(nanmean(pop_PPC2Pul_vid))));
mnPul2PPC = flipud(rot90(squeeze(nanmean(pop_Pul2PPC_vid))));

figure;
subplot(1,2,1)
% PPC 2 Pulvinar
imagesc(stepCen,1:numel(logF),mnPPC2Pul);
set(gca,'TickDir','out','YTick',tickLoc,'YDir','normal')
ax = gca;
ax.YTickLabel = {'0.5','1','2','4','8','16','32','64','128'};
xlabel('Time to saccade (s)')
ylabel('Frequency')
colormap(awesomemap2)
caxis([0 0.07])
ylim([200 800])
title('PPC to Pul video')
h = colorbar;
ylabel(h,'Granger causality')

% Pulvinar 2 PPC
subplot(1,2,2)
% PPC 2 Pulvinar
imagesc(stepCen,1:numel(logF),mnPul2PPC);
set(gca,'TickDir','out','YTick',tickLoc,'YDir','normal')
ax = gca;
ax.YTickLabel = {'0.5','1','2','4','8','16','32','64','128'};
xlabel('Time to saccade (s)')
ylabel('Frequency')
colormap(awesomemap2)
caxis([0 0.1])
ylim([200 800])
title('Pul to PPC video')
h = colorbar;
ylabel(h,'Granger causality')


% do stats
refSamp  = find(stepCen==-10);
testSamp = find(stepCen==0);
x = mean(pop_PPC2Pul(:,refSamp,al),3);
y = mean(pop_PPC2Pul(:,testSamp,al),3);
[h,p,ci,stats] = ttest(x,y);

% for non-parametric test
h  = kstest(x)
v  = vertcat(x,y);
gp = vertcat(ones(size(x)),ones(size(y))*2);
p  = kruskalwallis(v,gp)
[p,h,stats]  = signtest(x,y)

% plot PPC to Pul theta timecourse
al = find(logF>3.5 & logF < 4.5);
mn = mean(mean(pop_Pul2PPC(:,:,al),3));
er = std(mean(pop_Pul2PPC(:,:,al),3))/sqrt(size(pop_Pul2PPC,1));
% errorbar(stepCen,mn,er,'k'); hold on
% plot(stepCen,mn,'ok')
plot(stepCen,mn,'k'); hold on
plot(stepCen,mn+er,'k--');
plot(stepCen,mn-er,'k--');
set(gca,'TickDir','out')

% do stats
refSamp  = find(stepCen==-10);
testSamp = find(stepCen==2);
x = mean(pop_PPC2Pul(:,refSamp,al),3);
y = mean(pop_PPC2Pul(:,testSamp,al),3);
[h,p,ci,stats] = ttest(x,y);

% for non-parametric test
h  = kstest(x)
v  = vertcat(x,y);
gp = vertcat(ones(size(x)),ones(size(y))*2);
p  = kruskalwallis(v,gp)
p  = signtest(x,y)

refSamp = find(stepCen==-10);
bs      = mnPul2PPC(:,refSamp);
bsMat   = repmat(bs,1,numel(stepCen));
subMat  = (mnPul2PPC-bsMat);
imagesc(stepCen,1:numel(logF),subMat);
set(gca,'TickDir','out','YTick',tickLoc,'YDir','normal')
ax = gca;
ax.YTickLabel = {'0.5','1','2','4','8','16','32','64','128'};
xlabel('Time to saccade (s)')
ylabel('Frequency')
colormap(bry_map)
caxis([-0.03 0.03])
ylim([200 800])






