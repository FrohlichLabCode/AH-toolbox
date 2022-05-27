clear

addpath('E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/');
skipRec = 1;
% set linear or log plotting scale
[foi, tickLoc, tickLabel] = getFoiLabel(2, 32, 100, 2); % (lowFreq, highFreq, numFreqs, linORlog)
window   = 1*1024; %about 1sec
noverlap = round(15*window/16);
lfpFs = 1000;

animalCodes = {'0171'};

for iAnimal = 1:numel(animalCodes)
    animalCode = animalCodes{iAnimal};
    OptoDir       = ['E:/FerretData/' animalCode '/OptoMask/'];
    fileInfo      = dir([OptoDir 'randomIPI_*.mat']);  %_Level6* detect files to load/convert  '_LateralVideo*'
    for irec = 1:numel(fileInfo)    
        recName = fileInfo(irec).name;   %recName = 'randomIPI_15Hz.mat';
        saveName = [recName(1:end-4) '_fft'];
        if exist([OptoDir 'spec/' saveName '.mat']) && skipRec == 1; continue; end
        splitName = regexp(recName, '[_.]', 'split'); % split on multiple deliminators
        freqName = splitName{2};
        dutyName = splitName{3};        
        targetFreq = freqName(1:end-2);
        TTLarray = is_load([OptoDir recName],'randomTTLarray');    
        
        for iChn = 1:size(TTLarray)
            [s(iChn,:,:),f,t] = spectrogram(TTLarray(iChn,:),window,noverlap,foi,lfpFs);            
        end
        pow = abs(s).^2;
        
        AH_mkdir([OptoDir 'spec/']);
        save([OptoDir 'spec/' saveName '.mat'],'pow','foi','t');
    end
        
    
stimFs = [5.2, 15];
stimNames = {'Theta','Alpha'};
dutyCycles = 50;%[10:10:40];
% plot IPI spectrogram by dutyCycle
fig = AH_figure(2,2,'IPI spec');
for iDuty = 1:numel(dutyCycles)
    dutyCycle = dutyCycles(iDuty);
    
for iFreq = 1:numel(stimFs)
    stimF = stimFs(iFreq);    
    subplot(2,2,iFreq)
    recName = ['regularIPI_duty' num2str(dutyCycle)];
    [regularTTLarray{1}, regularTTLarray{2}] = is_load([OptoDir recName '.mat'],'thetaTTLarray','alphaTTLarray');
    [s{iFreq},f,t] = spectrogram(regularTTLarray{iFreq},window,noverlap,foi,lfpFs);
    imagesc(t,1:numel(foi),pow2db(abs(s{iFreq}.^2))); % mean of channels
    title(['regularIPI ' num2str(stimF) 'Hz duty' num2str(dutyCycle)]);
    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel);
    colorbar;colormap(jet);xlabel('Time [s]');ylabel('Freq [Hz]');caxis([-40,40]);
    
    subplot(2,2,2+iFreq)
    recName = ['randomIPI_' num2str(stimF) 'Hz_duty' num2str(dutyCycle)];
    pow = is_load([OptoDir 'spec/' recName '_fft.mat'],'pow');  
    imagesc(t,1:numel(foi),squeeze(pow2db(mean(pow,1)))); % mean of channels
    title(['randomIPI ' num2str(stimF) 'Hz duty' num2str(dutyCycle)]);
    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel);
    colorbar;colormap(jet);     
    xlabel('Time [s]');%caxis([0,35]);
    ylabel('Freq [Hz]');       
end
saveName = ['IPISpec_duty' num2str(dutyCycle) '_fft'];
savefig(fig, [OptoDir 'spec/' saveName '.fig'],'compact');
saveas(fig, [OptoDir 'spec/' saveName '.png']);
end
end