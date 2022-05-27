function region_spec_by_trial(skipRec, linORlog, lowFreq, highFreq, numFreqs, lfpFs,...
    evtTimes,twins,baseTwins, condNames, condID, regionNames, ...
    regionLFP, regionChn, sessionName,GroupAnalysisDir)

% normally use regionPair_FunConn is sufficient (for LateralVideo and
% PulvOpto) This V2 version allows different time windows to be used in
% different conditions, i.e. in CSRTT_FunConn, we need to align processing
% time with different time points: initiation, stimOnset, touch, etc.


doPlot = 1;
numConds   = numel(condID);
numRegions = numel(regionNames);

% Define frequencies of interest. Linear spacing for Phase slope index, and spacing for all other methods.
if     linORlog == 1
        foi      = linspace(lowFreq,highFreq,numFreqs); % linear spacing
elseif linORlog == 2
        foi      = logspace(log10(lowFreq),log10(highFreq),numFreqs); % log spacing
end    
GroupSessionDir = [GroupAnalysisDir 'sessions/'];
for iRegion = 1:numel(regionNames)
    regionName = regionNames{iRegion};
    
    if ~exist(join(GroupSessionDir),'dir'); mkdir(join(GroupSessionDir));
    elseif length(dir([GroupSessionDir 'sessionSpec_' sessionName '*.fig'])) >= length(condID) %each condition generate a fig file
        fprintf('Record %s already analyzed \n',sessionName'); 
        if skipRec == 1; continue; end  % to skip already analyzed records
    end
    fprintf('\nWorking on session %s \n',sessionName');     
    
    regionxLFP = regionLFP{iRegion};
    numChnReg1 = numel(regionChn{iRegion});

    for iCond = 1:numConds
        %if length(dir([saveAnalysisDir '*.fig'])) >= iCond;continue; end % skip already processed condition
        
        evtTime = evtTimes{condID(iCond)};     
        twin    = twins{condID(iCond)};
        baseTwin = baseTwins{condID(iCond)};
       
        for iReg1Chn = 1:numChnReg1
            xser = regionxLFP(iReg1Chn,:);
            fprintf('\nWorking on channel %d/%d, OuterCounter=%d================================ \n',...
                iReg1Chn, numChnReg1);

            clear funcCon % to clear out memory
            funcCon   = spec_by_trial(xser,regionName,condNames{condID(iCond)},...
                lfpFs,evtTime,twin,baseTwin,regionChn{iRegion}(iReg1Chn),...
                GroupAnalysisDir, linORlog, lowFreq, highFreq, numFreqs);

%                 funcCon = is_functionalConnectivity_V2(xser,yser,regionXname,regionYname,condNames{iCond},...
%                     lfpFs,evtTime,twin,baseTwin,regionChn{regionXind}(iReg1Chn),regionChn{regionYind}(iReg2Chn),...
%                     saveAnalysisDir, linORlog, lowFreq, highFreq, numFreqs);
            close all

            % for each struct matrix, it's frequency by time
            Spec.(regionName)(iReg1Chn,:,:,:) = funcCon.xspec;
            Spec.([regionName '_normed'])(iReg1Chn,:,:,:) = funcCon.xspecNormed;              
            
        end

    %% Eliminate nan sessions
    % The results from some sessions have nan elements. Those sessions need to
    % be eliminated. Note that simple use of "nanmean" does not work at this stage, but will be used later (EN).

%     for i=1:size(Spec.(regionName), 1)
%         for j=1:size(Spec.(regionName), 2)
%             temp1 = squeeze(Spec.(regionName)(i, j, :, :));
%             temp2 = isnan(temp1);
%             if any(temp2(:))
%                 fprintf('\n[i, j]= [%d, %d], NaN in xSpec   ', i, j);
%                 Spec.(regionName)(i, j, :, :) = nan;
%             end
%             temp1 = squeeze(Spec.([regionName '_normed'])(i, j, :, :));
%             temp2 = isnan(temp1);
%             if any(temp2(:))
%                 fprintf('\n[i, j]= [%d, %d], NaN in xSpec   ', i, j);
%                 Spec.([regionName '_normed'])(i, j, :, :) = nan;
%             end
% 
%             clear temp1 temp2
%         end
%     end

    %% Compute the mean values across channel pairs
    tvec = funcCon.tvec;
    
    % Compute ticks for plotting
    


    % calculating median
    TrialSpec.(regionName) = squeeze(nanmedian(Spec.(regionName),1)); % change mean to nanmedian over channels
    TrialSpec.([regionName '_normed']) = squeeze(nanmedian(Spec.([regionName '_normed']),1));
%     SessionSpec.(regionName) = squeeze(nanmedian(nanmedian(Spec.(regionName),1),2)); 
%     SessionSpec.([regionName '_normed']) = squeeze(nanmedian(nanmedian(Spec.([regionName '_normed']),1),2));
    
    
%     save([GroupSessionDir 'chnTrialSpec_' sessionName '_' regionName '_' condNames{condID(iCond)} '.mat'],'tvec','foi','fois','tickLoc','tickLabel','Spec','-v7.3');
    save([GroupSessionDir 'TrialSpec_' sessionName '_' regionName '_' condNames{condID(iCond)} '.mat'],'tvec','foi','fois','tickLoc','tickLabel','TrialSpec','baseTwins','-v7.3');
%     save([GroupSessionDir 'SessionSpec_' sessionName '_' regionName '_' condNames{condID(iCond)} '.mat'],'tvec','foi','fois','tickLoc','tickLabel','SessionSpec','baseTwins','-v7.3');
    
    fprintf('\nDone saving median ============================================\n')
    clear Spec TrialSpec SessionSpec
    
    end
end
    
    % to save memory space and also avoid using old data
    clear Spec TrialSpec


    
%% plotting
    if doPlot == 1
        
        %% first plot based on median        
        
% Plot
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 240*2 180*numRegions]); %(x,y,width,height) screensize(3)-100 480 x 540

for iCond = 1:numConds
    for iRegion = 1:numRegions
        regionName = regionNames{iRegion};
        load([GroupSessionDir 'SessionSpec_' sessionName '_' regionName '_' condNames{condID(iCond)} '.mat'])
        % plot power spectrum for signal x
        subplot(numRegions,2,iRegion)
        imagesc(tvec,1:numel(foi),pow2db(SessionSpec.(regionName)));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]'); % title('Signal X power')
    %    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        %caxis([22 55]);
        cl = colorbar('northoutside'); ylabel(cl,['Power [dB]: ' regionName],'FontSize',12)

        % plot power spectrum for signal y
        subplot(numRegions,2,iRegion+1)
        imagesc(tvec,1:numel(foi),pow2db(SessionSpec.([regionName '_normed'])));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Signal Y power')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        %caxis([22 55]);
        cl = colorbar('northoutside'); ylabel(cl,['Normed power: ' regionName],'FontSize',12)
    end
end
colormap(jet)

savefig(fig, [GroupSessionDir 'sessionSpec_' sessionName '_' condNames{condID(iCond)} '.fig'],'compact');
saveas(fig, [GroupSessionDir 'sessionSpec_' sessionName '_' condNames{condID(iCond)} '.png']);
close all
clear SessionSpec
end
end