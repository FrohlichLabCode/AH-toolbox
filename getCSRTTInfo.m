function [alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level)

delayNames   = {'D4','D5','D6','Dall'};
delayTypes   = {4,5,6,'Dall'};
optoNames    = {'Theta', 'Alpha', 'ArTheta','ArAlpha','Sham','Oall'};
hitMissNames = {'Correct','Premature','Incorrect','Omission','noPremature'};
% Trial inforamtion
if level(1) == '6'
    alignNames = {'Init','Stim','Touch'};
%     if length(level)>1 && level(2) == 'a'
%         delayNames = {'D2'};
%         delayTypes = {2};
%     end
else % level 7 8 9 
    alignNames = {'Init','Stim','Touch','Opto'};
    optoNames  = {'Theta','Alpha','ArTheta','ArAlpha','Sham','Oall'};
end
end