numFreqs = 100;
lowFreq  = 0.5;
highFreq = 128;
% foi      = linspace(lowFreq,highFreq,numFreqs); % linear spacing
foi   = linspace(lowFreq,highFreq,numFreqs); % log spacing

%% plotting
if 1
    
avgXSpec = squeeze(nanmean(nanmean( xSpecAll ,1),2));
avgYSpec = squeeze(nanmean(nanmean( ySpecAll ,1),2));
avgPLV = squeeze(nanmean(nanmean( plvAll ,1),2));
avgCoherencey = squeeze(nanmean(nanmean( coherencyAll ,1),2));
avgGC_XtoY = squeeze(nanmean(nanmean( GC_XtoY_All ,1),2));
avgGC_YtoX = squeeze(nanmean(nanmean( GC_YtoX_All ,1),2));
    
    
fois = 1:8:65;

% Compute ticks for plotting
    screensize = get( groot, 'Screensize' );
    fig = figure('Position',[10 50 screensize(3)-100 screensize(4)-150]);
    % Compute ticks for plotting
%    fois = [0.5 1 2 4 8 16 32 64 128];
    for fi = 1:numel(fois)
        [bi,bb] = sort(abs(foi-fois(fi)));
        tickLoc(fi) = bb(1);
    end
%    tickLabel = {'0.5','1','2','4','8','16','32','64','128'};
    tickLabel = {'2','10','18','26','34','42','50','58'};
    
    % do the same for phase slope index frequencies
%    fois = [0.6 1 2 4 8 16 32 64];
%     for fi = 1:numel(fois)
%         [bi,bb] = sort(abs(funcCon.psiFreq-fois(fi)));
%         psitickLoc(fi) = bb(1);
%     end
%     psitickLabel = {'0.6','1','2','4','8','16','32','64'};
    
    try
    % plot power spectrum for signal x
    subplot(2,4,1)
%     imagesc(tvec,1:numel(foi),pow2db(avgXSpec));
    imagesc(tvec,foi,pow2db(avgYSpec));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)'); % title('Signal X power')
%     set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    set(gca,'YDir','normal','TickDir','out')
    cl = colorbar('northoutside'); ylabel(cl,'Power (dB) signal X','FontSize',12)
    ylim([2,64]);
        catch
    end
    try
    % plot power spectrum for signal y
        subplot(2,4,2)
    imagesc(tvec,1:numel(foi),pow2db(avgYSpec));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('Signal Y power')
    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    cl = colorbar('northoutside'); ylabel(cl,'Power (dB) signal Y','FontSize',12)
    catch
    end
    try
    % plot phase locking value
    subplot(2,4,3)
    imagesc(tvec,1:numel(foi),avgPLV);
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('PLV')
    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    cl = colorbar('northoutside'); ylabel(cl,'Phase locking value','FontSize',12)
    catch
    end
    
    try
        % plot coherence
    subplot(2,4,4)
    imagesc(tvec,1:numel(foi),abs(avgCoherencey));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('Coherence')
    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    cl = colorbar('northoutside'); ylabel(cl,'Coherence','FontSize',12)
   
    catch
    end
    
    try
        % plot imaginary coherence
    subplot(2,4,5)
    imagesc(tvec,1:numel(foi),imag(avgCoherencey));
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('Imaginary coherence')
    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    %cl = colorbar('northoutside'); ylabel(cl,'Imaginary coherence (z-score)','FontSize',15)
    cl = colorbar('northoutside'); ylabel(cl,'Imag coherence (z)','FontSize',12)
    % plot phase slope index
    
    catch
    end
    
    try
        subplot(2,4,6)
    imagesc(tvec,1:numel(psiFreq),avgpsiNorm);
    xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('Phase slope index')
    set(gca,'YDir','normal','TickDir','out','YTick',psitickLoc,'YTickLabel',psitickLabel)
    %cl = colorbar('northoutside'); ylabel(cl,'Phase slope index (z-score)','FontSize',15)
    cl = colorbar('northoutside'); ylabel(cl,'Phase slope index (z)','FontSize',12)
    % plot granger causality X to Y
    catch
    end
    try
        subplot(2,4,7)
        imagesc(tvecGC,1:numel(foi),real(avgGC_XtoY));
        xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: X to Y')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        %cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: X to Y','FontSize',15)
        cl = colorbar('northoutside'); ylabel(cl,'GC: X to Y','FontSize',12)
        caxis([0 0.3]); ylim([1 90]); % 90+ Hz has saturated values
        % plot granger causality Y to X
        subplot(2,4,8)
        imagesc(tvecGC,1:numel(foi),real(avgGC_YtoX));
        xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: Y to X')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        %cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: Y to X','FontSize',15)
        cl = colorbar('northoutside'); ylabel(cl,'GC: Y to X','FontSize',12)
        caxis([0 0.3]); ylim([1 90]) % 90+ Hz has saturated values
    catch
    end
    colormap(jet)
    
    savefig(fig, 'avgFuncCon.fig','compact');
    saveas(fig, 'avgFuncCon.png');
end
x = toc;
fprintf('time required =%f sec\n', x);