

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
segLength = 2;
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
        % cz_detectSaccadesEphys(recPath,animalCode);
        load([recPath 'saccades'])
        
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
        
        %         % compute ERP and subtract on a trial-by-trial basis
        %         erpWin = round([-5 15]*newFs);
        %         ierp   = nan(numel(picOnset),diff(erpWin)+1);
        %         jerp   = nan(numel(picOnset),diff(erpWin)+1);
        %         for iev = 1:numel(picOnset)
        %             samps = round(picOnset(iev)*newFs)+erpWin(1):round(picOnset(iev)*newFs)+erpWin(2);
        %             ierp(iev,:) = idat(samps);
        %             jerp(iev,:) = jdat(samps);
        %         end
        
        
        
        numSac = [];
        for iev = 1:numel(vidOnset)
            numSac(iev)    = sum(sacTime>vidOnset(iev)&sacTime<vidOnset(iev)+10);
        end
        vidSac = numSac;
        
        numSac = [];
        for iev = 1:numel(picOnset)
            numSac(iev)    = sum(sacTime>picOnset(iev)&sacTime<picOnset(iev)+10);
        end
        picSac = numSac;
        
        segLength = 1;
        recLength = numel(idat)/newFs;
        nobs      = round(segLength*newFs);
        
        % Images
        
        upSamps = [];
        downSamps = [];
        for ipic = 1:numel(picOnset)
            samps = round(picOnset(ipic)*newFs):round(picOnset(ipic)*newFs)+(10*newFs)-1;
            if picSac(ipic) >= 3
                upSamps = [upSamps samps];
            else
                downSamps = [downSamps samps];
            end
        end
        
        % do down first
        imat = reshape(idat(downSamps),[nobs round(numel(downSamps)/newFs)]);
        jmat = reshape(jdat(downSamps),[nobs round(numel(downSamps)/newFs)]);
        
        clear X
        X(1,:,:)  = imat;
        X(2,:,:)  = jmat;
        numSeg    = size(imat,2);
        
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
        
        var_info(info,true); % report results (and bail out on error)
        
        F = autocov_to_pwcgc(G);
        assert(~isbad(F,false),'GC calculation failed');
        
        pval = mvgc_pval(F,morder,nobs,numSeg,1,1,nvars-2,tstat); % take careful note of arguments!
        sig  = significance(pval,alpha,mhtc);
        
        f = autocov_to_spwcgc(G,fres);
        assert(~isbad(f,false),'spectral GC calculation failed');
        freqRes        = size(f,3)-1;
        freqs          = linspace(0,newFs/2,freqRes+1)';
        
        PPC2Pul_down_pic(count,:) = interp1(freqs,squeeze(f(2,1,:)),logF,'spline');
        Pul2PPC_down_pic(count,:) = interp1(freqs,squeeze(f(1,2,:)),logF,'spline');
        
        
        if ~isempty(upSamps)
            % now do up
            imat = reshape(idat(upSamps),[nobs round(numel(upSamps)/newFs)]);
            jmat = reshape(jdat(upSamps),[nobs round(numel(upSamps)/newFs)]);
            
            clear X
            X(1,:,:)  = imat;
            X(2,:,:)  = jmat;
            numSeg    = size(imat,2);
            
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
            
            [A,SIG] = tsdata_to_var(X,morder,regmode);
            assert(~isbad(A),'VAR estimation failed');
            
            [G,info] = var_to_autocov(A,SIG,acmaxlags);
            
            var_info(info,true); % report results (and bail out on error)
            
            F = autocov_to_pwcgc(G);
            assert(~isbad(F,false),'GC calculation failed');
            
            pval = mvgc_pval(F,morder,nobs,numSeg,1,1,nvars-2,tstat); % take careful note of arguments!
            sig  = significance(pval,alpha,mhtc);
            
            f = autocov_to_spwcgc(G,fres);
            assert(~isbad(f,false),'spectral GC calculation failed');
            freqRes        = size(f,3)-1;
            freqs          = linspace(0,newFs/2,freqRes+1)';
            
            PPC2Pul_up_pic(count,:) = interp1(freqs,squeeze(f(2,1,:)),logF,'spline');
            Pul2PPC_up_pic(count,:) = interp1(freqs,squeeze(f(1,2,:)),logF,'spline');
        end
        
        % Videos
        
        upSamps = [];
        downSamps = [];
        for ivid = 1:numel(vidOnset)
            samps = round(vidOnset(ivid)*newFs):round(vidOnset(ivid)*newFs)+(10*newFs)-1;
            if vidSac(ivid) >= 3
                upSamps = [upSamps samps];
            else
                downSamps = [downSamps samps];
            end
        end
        
        % do down first
        imat = reshape(idat(downSamps),[nobs round(numel(downSamps)/newFs)]);
        jmat = reshape(jdat(downSamps),[nobs round(numel(downSamps)/newFs)]);
        
        clear X
        X(1,:,:)  = imat;
        X(2,:,:)  = jmat;
        numSeg    = size(imat,2);
        
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
        
        var_info(info,true); % report results (and bail out on error)
        
        F = autocov_to_pwcgc(G);
        assert(~isbad(F,false),'GC calculation failed');
        
        pval = mvgc_pval(F,morder,nobs,numSeg,1,1,nvars-2,tstat); % take careful note of arguments!
        sig  = significance(pval,alpha,mhtc);
        
        f = autocov_to_spwcgc(G,fres);
        assert(~isbad(f,false),'spectral GC calculation failed');
        freqRes        = size(f,3)-1;
        freqs          = linspace(0,newFs/2,freqRes+1)';
        
        PPC2Pul_down_vid(count,:) = interp1(freqs,squeeze(f(2,1,:)),logF,'spline');
        Pul2PPC_down_vid(count,:) = interp1(freqs,squeeze(f(1,2,:)),logF,'spline');
        
        
        if ~isempty(upSamps)
            % now do up
            imat = reshape(idat(upSamps),[nobs round(numel(upSamps)/newFs)]);
            jmat = reshape(jdat(upSamps),[nobs round(numel(upSamps)/newFs)]);
            
            clear X
            X(1,:,:)  = imat;
            X(2,:,:)  = jmat;
            numSeg    = size(imat,2);
            
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
            
            [A,SIG] = tsdata_to_var(X,morder,regmode);
            assert(~isbad(A),'VAR estimation failed');
            
            [G,info] = var_to_autocov(A,SIG,acmaxlags);
            
            var_info(info,true); % report results (and bail out on error)
            
            F = autocov_to_pwcgc(G);
            assert(~isbad(F,false),'GC calculation failed');
            
            pval = mvgc_pval(F,morder,nobs,numSeg,1,1,nvars-2,tstat); % take careful note of arguments!
            sig  = significance(pval,alpha,mhtc);
            
            f = autocov_to_spwcgc(G,fres);
            assert(~isbad(f,false),'spectral GC calculation failed');
            freqRes        = size(f,3)-1;
            freqs          = linspace(0,newFs/2,freqRes+1)';
            
            PPC2Pul_up_vid(count,:) = interp1(freqs,squeeze(f(2,1,:)),logF,'spline');
            Pul2PPC_up_vid(count,:) = interp1(freqs,squeeze(f(1,2,:)),logF,'spline');
        end
    end
end

PPC2Pul_down_pic(6,:) = nan;
PPC2Pul_up_pic(6,:) = nan;
PPC2Pul_down_vid(6,:) = nan;
PPC2Pul_up_vid(6,:) = nan;
Pul2PPC_down_pic(6,:) = nan;
Pul2PPC_up_pic(6,:) = nan;
Pul2PPC_down_vid(6,:) = nan;
Pul2PPC_up_vid(6,:) = nan;
PPC2Pul_down_pic(abs(PPC2Pul_down_pic)==0) = nan;
PPC2Pul_up_pic(abs(PPC2Pul_up_pic)==0) = nan;
PPC2Pul_down_vid(abs(PPC2Pul_down_vid)==0) = nan;
PPC2Pul_up_vid(abs(PPC2Pul_up_vid)==0) = nan;
Pul2PPC_down_pic(abs(Pul2PPC_down_pic)==0) = nan;
Pul2PPC_up_pic(abs(Pul2PPC_up_pic)==0) = nan;
Pul2PPC_down_vid(abs(Pul2PPC_down_vid)==0) = nan;
Pul2PPC_up_vid(abs(Pul2PPC_up_vid)==0) = nan;

subplot(2,2,1)
semilogx(logF,nanmean(PPC2Pul_down_pic),'b'); hold on
semilogx(logF,nanmean(PPC2Pul_up_pic),'r')
xlim([foi(1) foi(end)])
legend('< 4 saccades','>= 4 saccades')
xlabel('Frequency')
ylabel('Granger causality')
title('PPC->Pul Images')

subplot(2,2,2)
semilogx(logF,nanmean(PPC2Pul_down_vid),'b'); hold on
semilogx(logF,nanmean(PPC2Pul_up_vid),'r')
xlim([foi(1) foi(end)])
legend('< 4 saccades','>= 4 saccades')
xlabel('Frequency')
ylabel('Granger causality')
title('PPC->Pul Videos')

subplot(2,2,3)
semilogx(logF,nanmean(Pul2PPC_down_pic),'b'); hold on
semilogx(logF,nanmean(Pul2PPC_up_pic),'r')
xlim([foi(1) foi(end)])
legend('< 4 saccades','>= 4 saccades')
xlabel('Frequency')
ylabel('Granger causality')
title('Pul->PPC Images')

subplot(2,2,4)
semilogx(logF,nanmean(Pul2PPC_down_vid),'b'); hold on
semilogx(logF,nanmean(Pul2PPC_up_vid),'r')
xlim([foi(1) foi(end)])
legend('< 4 saccades','>= 4 saccades')
xlabel('Frequency')
ylabel('Granger causality')
title('Pul->PPC Videos')

apeak = 593;
tpeak = 366;


subplot(2,2,1)
d(:,1) = PPC2Pul_down_pic(:,apeak);
d(:,2) = PPC2Pul_up_pic(:,apeak);
boxplot(d)
[h,p,ci,stats] = ttest(d(:,1),d(:,2));
set(gca,'XTickLabel',{'< 4 sac','>=4 sac'})
ylabel('Granger causality (alpha)')
title(['Images - ttest p val = ' num2str(p)])

subplot(2,2,2)
d(:,1) = PPC2Pul_down_vid(:,apeak);
d(:,2) = PPC2Pul_up_vid(:,apeak);
boxplot(d)
[h,p,ci,stats] = ttest(d(:,1),d(:,2));
set(gca,'XTickLabel',{'< 4 sac','>=4 sac'})
ylabel('Granger causality (alpha)')
title(['Video - ttest p val = ' num2str(p)])

subplot(2,2,3)
d(:,1) = Pul2PPC_down_pic(:,tpeak);
d(:,2) = Pul2PPC_up_pic(:,tpeak);
boxplot(d)
[h,p,~,stats] = ttest(d(:,1),d(:,2));
set(gca,'XTickLabel',{'< 4 sac','>=4 sac'})
ylabel('Granger causality (theta)')
title(['Images - ttest p val = ' num2str(p)])

subplot(2,2,4)
d(:,1) = Pul2PPC_down_vid(:,tpeak);
d(:,2) = Pul2PPC_up_vid(:,tpeak);
boxplot(d)
[h,p,ci,stats] = ttest(d(:,1),d(:,2));
set(gca,'XTickLabel',{'< 4 sac','>=4 sac'})
ylabel('Granger causality (theta)')
title(['Video - ttest p val = ' num2str(p)])



fois = [0.5 1 2 4 8 16 32 64 128];
for fi = 1:numel(fois)
    [bi,bb] = sort(abs(logF-fois(fi)));
    tickLoc(fi) = bb(1);
end

mnPPC2Pul = flipud(rot90(squeeze(nanmean(pop_PPC2Pul_pic))));
mnPul2PPC = flipud(rot90(squeeze(nanmean(pop_Pul2PPC_pic))));

figure;
subplot(2,2,1)
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
subplot(2,2,2)
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

subplot(2,2,3)
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

% Pulvinar 2 PPC
subplot(2,2,4)
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



% plot PPC to Pul alpha timecourse
al = find(logF>13 & logF < 17);
mn = mean(mean(pop_PPC2Pul(:,:,al),3));
er = std(mean(pop_PPC2Pul(:,:,al),3))/sqrt(size(pop_PPC2Pul,1));
% errorbar(stepCen,mn,er,'k'); hold on
% plot(stepCen,mn,'ok')
plot(stepCen,mn,'k'); hold on
plot(stepCen,mn+er,'k--');
plot(stepCen,mn-er,'k--');
set(gca,'TickDir','out')

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






