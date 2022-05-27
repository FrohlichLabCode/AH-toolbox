function [validChn, reorderedChn] = keepChn(recName)
% This function selects valid channels and assign channel numbers to the regions of interst.
% Inputs:  recName      - eg. "0168_Opto_005_20180710"
% Outputs: validChn     - eg. {[1,16],[19,31],[35,40]} a cell array of regions, 
%                           each region contains a list of sorted channel index 
%                           corresponding to the row # of the lfpMat matrix.
%          reorderedChn - eg. {[1,2],[3,4],[5,6]} a cell array of regions, 
%                           each region contrains a list of new chanel index, 
%                           up to total number of valid channels.
%
% To assign valid channels, one can either put in channel # directly or use the
% chanOrder list to get a better sense of the spatial location of the channels
% 
% Use: [lfp.validChn,~] = keepChn(recName);
%
% AH 6/30/2018 created
% AH 3/26/2021 added check if pul channel is annatomically valid for 0171 0179 0180 0181

splitName   = strsplit(recName,'_'); % eg. output {'0168'}    {'Opto'}    {'005'}    {'20180710'}
animalCode  = splitName{1};
sessionType = splitName{2};
sessionID   = splitName{3};
date        = splitName{4}(1:8);

chanOrder16 = [5 6 7 8 9 10 11 12 ; 4 3 2 1 16 15 14 13]; % medial to lateral
chanOrder32 = [24 25 26 27 28 29 30 31 0 1 2 3 4 5 6 7 ; 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8] +1; % medial to lateral (word facing anterior, first row = posterior row.

switch animalCode
    case {'0201'}
        pmcChn = [1:16]; %201911
        claChn = [1:16] + 16;
        ppcChn = [1:16] + 32;
        if strcmp(date, '20211216')
            pmcChn = [1,4:6,11:14,16];
            claChn = [1:4,14:16]+16;
            ppcChn = [1,3,5,6,12,13,15:16] + 32;
        end
        validChn = {pmcChn; claChn; ppcChn};
    case {'0187', '0188'}
        ppcChn = [1:16];
        vcChn  = [1:10,12:16] + 16;
        validChn = {ppcChn; vcChn};
    case {'0179'}
        fcChn  = [1:16]; %201911
        pulChn = [1:16] + 16;
        ppcChn = [1:16] + 32;
        vcChn  = [1:16] + 48;
        if strcmp(date, '20191108')
            pulChn = [1:11,14:16] + 16;
            fcChn  = [1:3,6:7,15:16]; 
            ppcChn = [1:3,12,15:16] + 32;
            vcChn  = [1:7,9:10,13:16] + 48;
        elseif strcmp(date, '20191122')
            pulChn = [1:16] + 16;
            fcChn  = [1:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;
        elseif strcmp(date, '20191127')
            pulChn = [1:2,4,10:16] + 16;
            fcChn  = [1:16]; 
            ppcChn = [1:2,5,8:16] + 32;
            vcChn  = [1:16] + 48;
        elseif strcmp(date, '20191128')
            pulChn = [1:2,4,9:16] + 16;
            fcChn  = [1,3:6,9:16]; 
            ppcChn = [1:2,4:9,11:16] + 32;
            vcChn  = [1:16] + 48;
        elseif strcmp(date, '20191129') && strcmp(sessionID, '18')
            pulChn = [1,4,10:16] + 16;
            fcChn  = [1:7,10:16]; 
            ppcChn = [1:3,5,11:16] + 32;
            vcChn  = [1:14,16] + 48;
        elseif strcmp(date, '20191129') &&  strcmp(sessionID, '19')
            pulChn = [1,4,9:12,14:16] + 16;
            fcChn  = [1:7,10:11,13:16]; 
            ppcChn = [1:5,10:16] + 32;
            vcChn  = [1,3:12,14,16] + 48;
        elseif strcmp(date, '20191203')
            pulChn = [1,4,9:12,15:16] + 16;
            fcChn  = [1:7,9:11,13:16]; 
            ppcChn = [1:3,5,10:16] + 32;
            vcChn  = [1,3:6,11:12,14,16] + 48;
        elseif strcmp(date, '20191204')
            pulChn = [1,4,9:12,15:16] + 16;
            fcChn  = [1:7,9:11,13:16]; 
            ppcChn = [1:3,5,10:16] + 32;
            vcChn  = [1,3:6,11:12,14,16] + 48;
        elseif strcmp(date, '20191205')
            pulChn = [1:16] + 16;
            fcChn  = [1:7,9:11,13:16]; 
            ppcChn = [1:4,6:16] + 32;
            vcChn  = [1,3:6,9,11:13,16] + 48;
        elseif strcmp(date, '20191206')
            pulChn = [1:16] + 16;
            fcChn  = [1,3:7,9:11,13:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1,3:6,9,11:13,16] + 48;
        elseif strcmp(date, '20191209')
            pulChn = [1,4:16] + 16;
            fcChn  = [1,3:7,9:11,13,16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1,3:16] + 48;
        elseif strcmp(date, '20191211')
            pulChn = [1:16] + 16;
            fcChn  = [1,3:7,9:11,13,16]; 
            ppcChn = [1:4,6:16] + 32;
            vcChn  = [1,5:16] + 48;
        elseif strcmp(date, '20191211') && strcmp(sessionID, '01')
            pulChn = [1,5:16] + 16;
            fcChn  = [1,4:7,9:11,13,16]; 
            ppcChn = [1:4,6:14,16] + 32;
            vcChn  = [1,3,5:16] + 48;
        elseif strcmp(date, '20191212')
            pulChn = [1:16] + 16; %noisy
            fcChn  = [1,4:7,9:11,13,16]; 
            ppcChn = [1:4,6:16] + 32;
            vcChn  = [1:5,9,11:16] + 48;
        elseif strcmp(date, '20191213')
            pulChn = [1:12,14:16] + 16; 
            fcChn  = [1,4:5,7,10,16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:2,4:5,7,9,11:16] + 48;
        elseif strcmp(date, '20191216') && strcmp(sessionID, '05')
            pulChn = [1:3,5:12,14:16] + 16; 
            fcChn  = [1,4,7:8,10,16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:2,4:5,7,9,11:16] + 48;
        elseif strcmp(date, '20191216') && strcmp(sessionID, '06')
            pulChn = [1:12,14:16] + 16; 
            fcChn  = [1:12,14:16];  %new PFC headstage
            ppcChn = [1:16] + 32;
            vcChn  = [1,5,7,9,11:13,16] + 48;
        elseif strcmp(date, '20191217') && strcmp(sessionID, '26')
            pulChn = [1:6,8:9,12,14:16] + 16; 
            fcChn  = [1:12,14:16];  
            ppcChn = [1,3:16] + 32;
            vcChn  = [1,4,9,15,16] + 48;
        elseif strcmp(date, '20191217') && strcmp(sessionID, '07')
            pulChn = [1:9,11:12,14:16] + 16; 
            fcChn  = [1:12,14:16];  
            ppcChn = [1,3:16] + 32;
            vcChn  = [1,5:6,9,11:13,16] + 48;
        elseif strcmp(date, '20191218')
            pulChn = [1:9,11:12,14:16] + 16; 
            fcChn  = [1:12,14:16];  
            ppcChn = [1:16] + 32;
            vcChn  = [1,3:16] + 48;
        elseif strcmp(date, '20191219') %both sessions
            pulChn = [1:9,12:16] + 16; 
            fcChn  = [1:12,14:16];  
            ppcChn = [1,3:16] + 32;
            vcChn  = [1,3:16] + 48;
        elseif strcmp(date, '20191221') && strcmp(sessionID, '03')
            pulChn = [1:9,12,14] + 16; 
            fcChn  = [1:16];  
            ppcChn = [1,3:16] + 32;
            vcChn  = [1,3:16] + 48;            
        elseif strcmp(date, '20191221') && strcmp(sessionID, '04')
            pulChn = [1:9,12:16] + 16; 
            fcChn  = [1:3,5:16];  
            ppcChn = [1:16] + 32;
            vcChn  = [1,3:16] + 48;              
        elseif strcmp(date, '20191222') && strcmp(sessionID, '10')
            pulChn = [1:9,12,14:16] + 16; 
            fcChn  = [1:11,14:16]; 
            ppcChn = [1:15] + 32;
            vcChn  = [1:9,11:16] + 48;
        elseif strcmp(date, '20191222') && strcmp(sessionID, '11')
            pulChn = [1:9,12:16] + 16; 
            fcChn  = [1:11,14:16];  
            ppcChn = [1:6,8:12,14:16] + 32;
            vcChn  = [1:16] + 48;
        elseif strcmp(date, '20191223')
            pulChn = [1:3,5:16] + 16; 
            fcChn  = [1:3,5:11]; 
            ppcChn = [1:6,8:13] + 32;
            vcChn  = [1:16] + 48;    
        elseif strcmp(date, '20191228')
            pulChn = [1:16] + 16; 
            fcChn  = [1:6,8:16];  %new PFC headstage
            ppcChn = [1:5,7:16] + 32;
            vcChn  = [1,3:16] + 48;
        elseif strcmp(date, '20191229')
            pulChn = [1:16] + 16; 
            fcChn  = [1:16];  %new PFC headstage
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;        
        elseif strcmp(date, '20200102')
            pulChn = [1:3,5:14] + 16; 
            fcChn  = [1:16];  %new PFC headstage
            ppcChn = [1:16] + 32;
            vcChn  = [1,3:6,8:16] + 48;    
        elseif strcmp(date, '20200103')
            pulChn = [1:3,5:14] + 16; 
            fcChn  = [2:16];  %new PFC headstage
            ppcChn = [1:16] + 32;
            vcChn  = [1,3:6,8:16] + 48; % -- above all filled
        elseif strcmp(date, '20200109') % both
            pulChn = [1:16] + 16; 
            fcChn  = [1:3,5:16];  %new PFC headstage
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48; 
        elseif strcmp(date, '20200110')
            pulChn = [1:7,9:16] + 16; 
            fcChn  = [1:12,14:16];  %new PFC headstage
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48; 
            
        elseif strcmp(date, '20200113') && strcmp(sessionID, '18')
            pulChn = [1:7,9:16] + 16; 
            fcChn  = [1:5,9:16];  %new PFC headstage
            ppcChn = [1:8,10:16] + 32;
            vcChn  = [1:16] + 48;    
        elseif strcmp(date, '20200113') 
            pulChn = [1:7,9:16] + 16; 
            fcChn  = [1:12,14:16];  %new PFC headstage
            ppcChn = [2:16] + 32;
            vcChn  = [1:16] + 48;
        elseif strcmp(date, '20200114')
            pulChn = [1:4,6,9:16] + 16; 
            fcChn  = [1:5,9:12,14:16];  %new PFC headstage
            ppcChn = [1:16] + 32;
            vcChn  = [1:6,8:16] + 48;
        elseif strcmp(date, '20200115') %both
            pulChn = [1:4,6,8:16] + 16; 
            fcChn  = [1:5,9:12,14:16];  
            ppcChn = [1:16] + 32;
            vcChn  = [1:6,8:16] + 48;
        elseif strcmp(date, '20200117') % both
            pulChn = [1:16] + 16; 
            fcChn  = [1:3,5,16];  
            ppcChn = [1:15] + 32;
            vcChn  = [1:6,8:16] + 48;
        elseif strcmp(date, '20200121')
            pulChn = [1:10,12:16] + 16; 
            fcChn  = [1:3,5,16]; 
            ppcChn = [1:7,9:16] + 32;
            vcChn  = [1:3,5:14,16] + 48;
        elseif strcmp(date, '20200122') % both
            pulChn = [1:16] + 16; 
            fcChn  = [1:3,5,16];  
            ppcChn = [1:16] + 32;
            vcChn  = [1:3,5:15] + 48;
        elseif strcmp(date, '20200123')
            pulChn = [1:16] + 16; 
            fcChn  = [1:16];  
            ppcChn = [1:4,6:16] + 32;
            vcChn  = [1:3,5:13] + 48;
        elseif strcmp(date, '20200124') 
            pulChn = [1:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [1:4,6:16] + 32;
            vcChn  = [1:13,16] + 48;            
        elseif strcmp(date, '20200127')
            pulChn = [1:16] + 16; 
            fcChn  = [1:16];  
            ppcChn = [1:16] + 32;
            vcChn  = [1:13,15:16] + 48;
        elseif strcmp(date, '20200128')
            pulChn = [1:16] + 16; 
            fcChn  = [1:16];  
            ppcChn = [1:5,7:16] + 32;
            vcChn  = [1:10,12:14,16] + 48; 
        elseif strcmp(date, '20200129')
            pulChn = [1:16] + 16; 
            fcChn  = [1:16];  %new PFC headstage
            ppcChn = [1:16] + 32;
            vcChn  = [1:11,13,14,16] + 48; % -- above all filled 1/30/2020
        elseif strcmp(date, '20200130') && strcmp(sessionID, '37')
            pulChn = [1:16] + 16; 
            fcChn  = [1:16];  %new PFC headstage
            ppcChn = [1:5,7:12,14:16] + 32;
            vcChn  = [1:2,5:14,16] + 48; 
        elseif strcmp(date, '20200130')
            pulChn = [1:16] + 16; 
            fcChn  = [1:16];  %new PFC headstage
            ppcChn = [1:5,7,10:15] + 32;
            vcChn  = [1:3,5:11,13,14,16] + 48;      
        elseif strcmp(date, '20200131')
            pulChn = [1:16] + 16; 
            fcChn  = [1:16];  %new PFC headstage
            ppcChn = [1:3,7:16] + 32;
            vcChn  = [1,5:14,16] + 48;  % -- above all filled 2/1/2020
        elseif strcmp(date, '20200203')
            pulChn = [1:16] + 16; 
            fcChn  = [1:16];  %new PFC headstage
            ppcChn = [1:3,5,7,9:16] + 32;
            vcChn  = [1:2,5:13,16] + 48;
        elseif strcmp(date, '20200204') % PFC noisy cable, not usable
            pulChn = [1:16] + 16; 
            fcChn  = [1:16];  %new PFC headstage
            ppcChn = [1:5,7,9:16] + 32;
            vcChn  = [1:2,4:12,16] + 48; 
        elseif strcmp(date, '20200206')
            pulChn = [1:16] + 16; 
            fcChn  = [1:4,7:14,16];  
            ppcChn = [1:2,4,7:12,14:16] + 32;
            vcChn  = [1:16] + 48;
        elseif strcmp(date, '20200207')
            pulChn = [1:3,5:9,12:16] + 16; 
            fcChn  = [1:10,12:16];  
            ppcChn = [1:2,4,6:7,9:10,12:16] + 32;
            vcChn  = [1:16] + 48;
        elseif strcmp(date, '20200214')
            pulChn = [1:5,7:9,12:16] + 16; 
            fcChn  = [1:7,9,13,15];  
            ppcChn = [1:5,7,9,12:16] + 32;
            vcChn  = [1:5,9:16] + 48;  % -- above all filled 2/17/2020
        elseif strcmp(date, '20200218')
            pulChn = [1:8,12:16] + 16; 
            fcChn  = [1:6,9,12:13,15:16];  
            ppcChn = [1:4,9,12:13,16] + 32;
            vcChn  = [1:4,8:16] + 48; 
        elseif strcmp(date, '20200219')
            pulChn = [1,3:5,7:9,12:16] + 16; 
            fcChn  = [1:6,9,13,15:16];  
            ppcChn = [1:4,9,12:16] + 32;
            vcChn  = [1:4,8:12,15:16] + 48; 
        elseif strcmp(date, '20200220')
            pulChn = [1,3:8,13:16] + 16; 
            fcChn  = [1,3:6,9,13,16];  
            ppcChn = [1:4,12:16] + 32;
            vcChn  = [1:4,9,11:12,16] + 48; 
        elseif strcmp(date, '20200225')
            pulChn = [1,3:9,12:16] + 16; 
            fcChn  = [1,3,9,13,16];  
            ppcChn = [1:4,12:16] + 32;
            vcChn  = [1:4,9,11:12,16] + 48; % -- above all filled 3/4/2020
        elseif strcmp(date, '20200305')
            pulChn = [1,7:10,12:14,16] + 16; 
            fcChn  = [1:16];  
            ppcChn = [1:4,12:14,16] + 32;
            vcChn  = [1:4,9,11:13,16] + 48;
        elseif strcmp(date, '20200306')
            pulChn = [1,7:8,10,12:14,16] + 16; 
            fcChn  = [1:9,12,14:16];  
            ppcChn = [1:4,12:14,16] + 32;
            vcChn  = [1:4,9,11:13,16] + 48; % -- above all filled 3/9/2020
              
        end
        %validChn = {fcChn; pulChn; ppcChn; vcChn};
        % Check if pul channel is annatomically valid (added 3/26/2021)
        validSUMask = getValidSUMask(animalCode,'LPl',pulChn-16);
        newpulChn = pulChn(validSUMask);
        validChn = {fcChn; newpulChn; ppcChn; vcChn};
%% 0180
    case {'0180'}
        fcChn  = [1:16]; %201911
        pulChn = [1:16] + 16;
        ppcChn = [1:16] + 32;
        vcChn  = [1:16] + 48;
        if strcmp(date, '20191025')
            fcChn  = [1:16]; 
            pulChn = [1:2,5:6,8:12,16] + 16;%line 20Hz
            ppcChn = [1:7,9:12,14:16] + 32;
            vcChn  = [1:5,13:16] + 48;
        elseif strcmp(date, '20191028')
            fcChn  = [1:16]; 
            pulChn = [1,4,6,8:12,16] + 16;%line 20Hz
            ppcChn = [1:9,11:12,14:16] + 32;
            vcChn  = [1:5,13:16] + 48;
        elseif strcmp(date, '20191029')
            fcChn  = [1:16]; 
            pulChn = [1,6,8:12,16] + 16;%line 20Hz
            ppcChn = [1:6,9,14:16] + 32;
            vcChn  = [1:5,13:16] + 48;
        elseif strcmp(date, '20191030')
            fcChn  = [1:16]; 
            pulChn = [1,6,8:12,16] + 16;%line 20Hz
            ppcChn = [1:2,4:9,12,14:16] + 32;
            vcChn  = [1:5,13:16] + 48;
        elseif strcmp(date, '20191105')
            fcChn  = [1,3,8,11:14,16]; 
            pulChn = [1,6:11,16] + 16;%line 20Hz
            ppcChn = [1:5,12,15:16] + 32;
            vcChn  = [1:3,5:6,10,13:16] + 48;
        elseif strcmp(date, '20191106')
            fcChn  = [1,3,8,11:14,16]; %noisy
            pulChn = [1:16] + 16;%noisy %line 20Hz
            ppcChn = [1:5,9,12,15:16] + 32;
            vcChn  = [1:5,10,13:16] + 48;            
        elseif strcmp(date, '20191107')
            fcChn  = [1,4,12,16]; %noisy
            pulChn = [1,4,6,10,12,16] + 16;%noisy %line 20Hz
            ppcChn = [1:2,4:5,9,12,15:16] + 32;
            vcChn  = [1:4,13,15:16] + 48;  
        elseif strcmp(date, '20191108')
            fcChn  = [1:3,14:16]; %noisy
            pulChn = [1:16] + 16;%noisy %line 20Hz
            ppcChn = [1:4,12,15:16] + 32;
            vcChn  = [1:16] + 48; 
        elseif strcmp(date, '20191111')
            fcChn  = [1:3,15:16]; %noisy
            pulChn = [1:5,10:11,15:16] + 16; %line 20Hz
            ppcChn = [1:3,15:16] + 32;
            vcChn  = [1:16] + 48; 
        elseif strcmp(date, '20191112')
            fcChn  = [1:3,14:16]; %noisy
            pulChn = [1:4,10:11,15:16] + 16; %noisy %line 20Hz
            ppcChn = [1:3,5,15:16] + 32;
            vcChn  = [1:16] + 48; %noisy
        elseif strcmp(date, '20191113')
            fcChn  = [1:3,16]; %noisy
            pulChn = [1:4,7,10:11,15:16] + 16; %noisy %line 20Hz
            ppcChn = [1:3,7,15:16] + 32;
            vcChn  = [1:16] + 48;            
        elseif strcmp(date, '20191114')
            fcChn  = [1:3,16]; %noisy
            pulChn = [1:5,10:11,15:16] + 16; %noisy %line 20Hz
            ppcChn = [1:2,7,15:16] + 32;
            vcChn  = [1:16] + 48; 
        elseif strcmp(date, '20191115') && strcmp(sessionID, '25')
            fcChn  = [1,3,16]; %noisy
            pulChn = [1:5,10:11,15:16] + 16; %noisy %line 20Hz
            ppcChn = [1:3,15:16] + 32;
            vcChn  = [1:16] + 48; 
        elseif strcmp(date, '20191115') && strcmp(sessionID, '26')
            fcChn  = [1:2,6:7,12,14,16]; %noisy
            pulChn = [1:5,10:11,15:16] + 16; %noisy %line 20Hz
            ppcChn = [1:2,4,15:16] + 32;
            vcChn  = [1:16] + 48; 
        elseif strcmp(date, '20191118')
            fcChn  = [1:2,5:7,12,14:16]; %noisy
            pulChn = [1:4,10:11,15:16] + 16; %noisy %line 20Hz
            ppcChn = [1:4,15:16] + 32;
            vcChn  = [1:16] + 48; 
        elseif strcmp(date, '20191119')
            fcChn  = [1:2,12,16]; %noisy
            pulChn = [1:4,10:11,15:16] + 16; %noisy %line 20Hz
            ppcChn = [1:4,15:16] + 32;
            vcChn  = [1:16] + 48; 
        elseif strcmp(date, '20191120')
            fcChn  = [1:2,16]; %noisy
            pulChn = [1:4,10:11,15:16] + 16; %noisy %line 20Hz
            ppcChn = [1:2,4,9,15:16] + 32;
            vcChn  = [1:16] + 48; 
        elseif strcmp(date, '20191121')
            fcChn  = [1:2,16]; %noisy
            pulChn = [1:4,10:11,15:16] + 16; %
            ppcChn = [1:2,4,15:16] + 32;
            vcChn  = [1:16] + 48; 
        elseif strcmp(date, '20191122')
            fcChn  = [1:2,16]; %noisy
            pulChn = [1:2,4:5,10:11,15:16] + 16; %
            ppcChn = [1:2,4:5,9,14:16] + 32;
            vcChn  = [1:16] + 48; 
        elseif strcmp(date, '20191125')
            fcChn  = [1:5,7:16]; %
            pulChn = [1:12,14:16] + 16; %
            ppcChn = [1:12,14,16] + 32;
            vcChn  = [1:10,12:16] + 48;
        elseif strcmp(date, '20191126')
            fcChn  = [1:16]; %
            pulChn = [1:8,10:12,14:16] + 16; %
            ppcChn = [1:9,11:16] + 32;
            vcChn  = [1:16] + 48;
        elseif strcmp(date, '20191127')
            fcChn  = [1:16]; %weird no signal
            pulChn = [1:5,10:12,14:16] + 16; %
            ppcChn = [1:2,5:6,11:16] + 32;
            vcChn  = [1:16] + 48;
        elseif strcmp(date, '20191128')
            fcChn  = [1:16]; %weird no signal
            pulChn = [1:2,4,9:12,14:16] + 16; %
            ppcChn = [1:2,4:6,8,11:16] + 32;
            vcChn  = [1:16] + 48;
        elseif strcmp(date, '20191129')
            fcChn  = [1:7,10:11,13:16]; %weird no signal
            pulChn = [1,4,9:16] + 16; %
            ppcChn = [1:5,9:16] + 32;
            vcChn  = [1:14,16] + 48;
        elseif strcmp(date, '20191203')
            fcChn  = [1:7,10:16]; %weird no signal
            pulChn = [1,4,9:12,15:16] + 16; %
            ppcChn = [1:3,5,10:16] + 32;
            vcChn  = [1,3:6,10:12,14,16] + 48;
        elseif strcmp(date, '20191204')
            fcChn  = [1,3:7,10:11,13:16]; 
            pulChn = [1,4,9:12,15:16] + 16; %
            ppcChn = [1:3,5,9:16] + 32;
            vcChn  = [1,3:6,9:12,14,16] + 48; 
        elseif strcmp(date, '20191206') % both
            fcChn  = [1,3:7,9:11,13:16]; 
            pulChn = [1,3:12,14:16] + 16; %
            ppcChn = [1:3,5:16] + 32;
            vcChn  = [1,3:6,10:13,16] + 48; 
        elseif strcmp(date, '20191209')
            fcChn  = [1,3:7,9:11,13:14,16]; 
            pulChn = [1,4:12,14:16] + 16; %
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;     
        elseif strcmp(date, '20191211')
            fcChn  = [1,3:5,7,10,13,16]; 
            pulChn = [1:12,14:16] + 16; %
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;
        elseif strcmp(date, '20191212')
            fcChn  = [1,4:5,7,10:11,13,16]; 
            pulChn = [1:12,14:16] + 16; %
            ppcChn = [1:16] + 32;
            vcChn  = [1:9,11:13,16] + 48;     
        elseif strcmp(date, '20191213')
            fcChn  = [1,4:5,7,10,13,16]; 
            pulChn = [1:12,14:16] + 16; %noisy
            ppcChn = [1:16] + 32;
            vcChn  = [1:2,4:5,9,11:13,16] + 48;   
        elseif strcmp(date, '20191213') && strcmp(sessionID, '03')
            fcChn  = [1,4:5,7,10,13,16]; 
            pulChn = [1:12,15:16] + 16;
            ppcChn = [1:10,12:16] + 32;
            vcChn  = [1:2,4:5,9,11:13,16] + 48;   
        elseif strcmp(date, '20191216') && strcmp(sessionID, '05')
            fcChn  = [1,4:5,7,10,16]; 
            pulChn = [1:12,14:16] + 16; %noisy
            ppcChn = [1:16] + 32;
            vcChn  = [1:2,4,7,9,11:16] + 48;               
        elseif strcmp(date, '20191216') && strcmp(sessionID, '04')
            fcChn  = [1:6,8:16]; 
            pulChn = [1:12,15:16] + 16; %noisy
            ppcChn = [1:16] + 32;
            vcChn  = [1,4:5,7,9,11:16] + 48; 
        elseif strcmp(date, '20191217')    
            fcChn  = [1:16]; 
            pulChn = [1:9,11:12,14:16] + 16; %noisy
            ppcChn = [1:16] + 32;
            vcChn  = [1,5,9,11,13,16] + 48;     
        elseif strcmp(date, '20191218')
            fcChn  = [1:16]; 
            pulChn = [1:7,9:12,15:16] + 16; %noisy
            ppcChn = [1:16] + 32;
            vcChn  = [1:12,14:16] + 48;   
        elseif strcmp(date, '20191219')
            fcChn  = [1:16]; 
            pulChn = [1:7,9:12,15:16] + 16; %noisy
            ppcChn = [1:16] + 32;
            vcChn  = [1:13,16] + 48;            
        elseif strcmp(date, '20191221')
            pulChn = [1:7,9,11:12,15:16] + 16; %noisy
            fcChn  = [1,3:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;            
        elseif strcmp(date, '20191221') && strcmp(sessionID, '06')
            pulChn = [1:7,9,11,13,15:16] + 16; %noisy
            fcChn  = [1,3:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;            
        elseif strcmp(date, '20191222')
            pulChn = [1:9,11:12,15:16] + 16; %noisy
            fcChn  = [1,3:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;
        elseif strcmp(date, '20191223')
            pulChn = [1:11,15:16] + 16; %noisy
            fcChn  = [1,3:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;            
        elseif strcmp(date, '20191228') %Lpl signal comes and goes
            pulChn = [1:11,15:16] + 16; 
            fcChn  = [1:12,14:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:9,11:12,16] + 48;             
        elseif strcmp(date, '20191229') %Lpl signal comes and goes
            pulChn = [1:12,15:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [1:3,5:16] + 32;
            vcChn  = [1:9,11:12,16] + 48;  
        elseif strcmp(date, '20200102') 
            pulChn = [1:12,14:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:6,8:16] + 48;             
        elseif strcmp(date, '20200103') 
            pulChn = [1:12,14:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:6,8:16] + 48;
        elseif strcmp(date, '20200108') 
            pulChn = [1:12,14:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:6,8:16] + 48;    
        elseif strcmp(date, '20200109') && strcmp(sessionID, '18')
            pulChn = [1:7,9:12,14:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:6,8:16] + 48;     
        elseif strcmp(date, '20200109') % both sessions
            pulChn = [1:7,9:12,14:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;    
        elseif strcmp(date, '20200110') 
            pulChn = [1:7,9:12,14:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;     
        elseif strcmp(date, '20200113') && strcmp(sessionID, '11')
            pulChn = [1:7,9:12,14:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;     
        elseif strcmp(date, '20200113')
            pulChn = [1:7,9:12,15:16] + 16; 
            fcChn  = [1:6,9:12,14:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;          
        elseif strcmp(date, '20200115') % both sessions
            pulChn = [1:4,6,8:12,14:16] + 16; 
            fcChn  = [1:5,9:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;        
        elseif strcmp(date, '20200116') 
            pulChn = [1:6,10:12,15:16] + 16; 
            fcChn  = [1:5,14:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:10,12:16] + 48;  
        elseif strcmp(date, '20200117') % both sessions
            pulChn = [1:11,16] + 16; 
            fcChn  = [1:3,5,16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:10,12:16] + 48;        
        elseif strcmp(date, '20200121') % both sessions
            pulChn = [1:12,15:16] + 16; 
            fcChn  = [1:3,5,16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;
        elseif strcmp(date, '20200122') && strcmp(sessionID, '35')
            pulChn = [1:7,9:11,15:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48; 
        elseif strcmp(date, '20200122')
            pulChn = [1:7,9:11,15:16] + 16; 
            fcChn  = [1:5,16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;  
        elseif strcmp(date, '20200123')
            pulChn = [1:7,9:11,15:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;    
        elseif strcmp(date, '20200124') % both sessions
            pulChn = [1:7,9:12,15:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [1:13,15:16] + 32;
            vcChn  = [1:14,16] + 48;
        elseif strcmp(date, '20200127') % both sessions
            pulChn = [1:7,9:12,15:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:11,13,15:16] + 48;
        elseif strcmp(date, '20200128')
            pulChn = [1:7,9:12,15:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [1:13,15:16] + 32;
            vcChn  = [1:13,15:16] + 48; % -- above all filled 1/30/2020
        elseif strcmp(date, '20200129') && strcmp(sessionID, '25')
            pulChn = [1,3:7,9:12,15:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [1:13,15:16] + 32;
            vcChn  = [1:13,16] + 48;       
        elseif strcmp(date, '20200129')
            pulChn = [1:11,15:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [1:13,15:16] + 32;
            vcChn  = [1:4,6:11,13,16] + 48; 
        elseif strcmp(date, '20200130') %%
            pulChn = [1:12,14:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [1:4:13,15:16] + 32;
            vcChn  = [1:4,6:11,13,16] + 48; 
        elseif strcmp(date, '20200131')
            pulChn = [1:12,15:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [1:13,15:16] + 32;
            vcChn  = [1:11,13:14,16] + 48; % -- above all filled 2/1/2020
        elseif strcmp(date, '20200203') % both
            pulChn = [1:7,8:11,15:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [1:4,6:7,9:13,15:16] + 32;
            vcChn  = [1:2,4:14,16] + 48; 
        elseif strcmp(date, '20200204')
            pulChn = [1:11,15:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [1:5,7:12,15:16] + 32;
            vcChn  = [1:2,4:13,16] + 48;
        elseif strcmp(date, '20200205')
            pulChn = [1:11,15:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [1:5,7:12,15:16] + 32;
            vcChn  = [1:2,4:13,16] + 48; 
        elseif strcmp(date, '20200206')
            pulChn = [1:10,15:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [2:5,7:12,15:16] + 32;
            vcChn  = [1:16] + 48;
        elseif strcmp(date, '20200207')
            pulChn = [1:10,15:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [1:3,6:7,9:10,12:13,15] + 32;
            vcChn  = [1:16] + 48; 
        elseif strcmp(date, '20200210')
            pulChn = [1:12,15:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [1:5,7,9:10,12:16] + 32;
            vcChn  = [1:16] + 48;     
        elseif strcmp(date, '20200213') % both
            pulChn = [1:9,12,15:16] + 16; 
            fcChn  = [1:12,14:16]; 
            ppcChn = [1:5,7,9,12:13,15:16] + 32;
            vcChn  = [1,3:16] + 48;  
        elseif strcmp(date, '20200214') % both
            pulChn = [1:5,7:9,11:12,15:16] + 16; 
            fcChn  = [1:6,9:10,12:13,15:16]; 
            ppcChn = [1:4,7,9,12:14,16] + 32;
            vcChn  = [1:6,8:16] + 48;  % -- above all filled 2/17/2020
        elseif strcmp(date, '20200217')
            pulChn = [1:12,14:16] + 16; 
            fcChn  = [1:10,12:13,15:16]; 
            ppcChn = [1:4,7,9,12,14:16] + 32;
            vcChn  = [1:4,6,8:13,16] + 48; 
        elseif strcmp(date, '20200218')
            pulChn = [1:12,14:16] + 16; 
            fcChn  = [1:6,9,12:13,15:16]; 
            ppcChn = [1:5,7,9,12:16] + 32;
            vcChn  = [1:4,6,8:13,16] + 48;
        elseif strcmp(date, '20200219')
            pulChn = [1,3:12,14:16] + 16; 
            fcChn  = [1:6,9,12:13,16]; 
            ppcChn = [1:5,7,9,12:16] + 32;
            vcChn  = [1:4,8:13,16] + 48;
        elseif strcmp(date, '20200220') && strcmp(sessionID, '01') % level9
            pulChn = [1,3:10,12,15:16] + 16; 
            fcChn  = [1,3,5:6,9,13,16]; 
            ppcChn = [1:4,7,9,12:14,16] + 32;
            vcChn  = [1:4,9,11:13,16] + 48; 
        elseif strcmp(date, '20200220') % both
            pulChn = [1,3:10,12,16] + 16; 
            fcChn  = [1,3,5:6,9,13,16]; 
            ppcChn = [1:4,7,12:14,16] + 32;
            vcChn  = [1:4,9,11:13,16] + 48;
        elseif strcmp(date, '20200224') && strcmp(sessionID, '34') %
            pulChn = [1,3:12,14,16] + 16; 
            fcChn  = [1,3,9,13,16]; 
            ppcChn = [1:5,7,9,12:16] + 32;
            vcChn  = [1:4,9,11:16] + 48;        
        elseif strcmp(date, '20200224') % both
            pulChn = [1,3:12,14,16] + 16; 
            fcChn  = [1,3,5:6,9,13,16]; 
            ppcChn = [1:5,7,9,12:16] + 32;
            vcChn  = [1:4,8:9,11:14,16] + 48;
        elseif strcmp(date, '20200225')
            pulChn = [1,3:12,14:16] + 16; 
            fcChn  = [1,3,9,13,16]; 
            ppcChn = [1:5,9,12:16] + 32;
            vcChn  = [1:4,9,11:13,16] + 48;
        elseif strcmp(date, '20200302')
            pulChn = [1,3:12,14:16] + 16; 
            fcChn  = [1,3,9,16]; 
            ppcChn = [1:5,7,12:16] + 32;
            vcChn  = [1:4,9,11:13,16] + 48;    
        elseif strcmp(date, '20200303') 
            pulChn = [1,3:10,12,16] + 16; 
            fcChn  = [1,9,16]; 
            ppcChn = [1:4,9,12:16] + 32; 
            vcChn  = [1:4,9,11:12,16] + 48; 
        elseif strcmp(date, '20200304') && strcmp(sessionID, '38')
            pulChn = [1,4:5,7:10,12,14,16] + 16; 
            fcChn  = [1:2,4:16]; 
            ppcChn = [1:4,7,12:13] + 32;
            vcChn  = [1:3,9,11:13,16] + 48;
        elseif strcmp(date, '20200304') 
            pulChn = [1,4:10,12,14,16] + 16; 
            fcChn  = [1:2,4:16]; 
            ppcChn = [1:4,12:13] + 32;
            vcChn  = [1:4,9,11:12,16] + 48;
        elseif strcmp(date, '20200305') && strcmp(sessionID, '03')
            pulChn = [1,7:12,14,16] + 16; 
            fcChn  = [1:2,4:9,11:16]; 
            ppcChn = [1:4,7,12:13] + 32;
            vcChn  = [1:4,9,11:12,16] + 48;  
        elseif strcmp(date, '20200305') 
            pulChn = [1,5,7:12,14,16] + 16; 
            fcChn  = [1:2,4:16]; 
            ppcChn = [1:4,12:13] + 32;
            vcChn  = [1:4,9,11:12,16] + 48;    
        elseif strcmp(date, '20200306') 
            pulChn = [1,7:10,12,14,16] + 16; 
            fcChn  = [1:2,4:9,11:16]; 
            ppcChn = [1:4,12:13] + 32;
            vcChn  = [1:4,9,11:12,16] + 48; % -- above all filled 3/9/2020
        end
        %validChn = {fcChn; pulChn; ppcChn; vcChn}; % end of 0180
        
        % Check if pul channel is annatomically valid (added 3/26/2021)
        validSUMask = getValidSUMask(animalCode,'LPl',pulChn-16);
        newpulChn = pulChn(validSUMask);
        validChn = {fcChn; newpulChn; ppcChn; vcChn};
%% 0181
    case '0181'
        fcChn  = [1:16]; %201911
        pulChn = [1:16] + 16;
        ppcChn = [1:16] + 32;
        vcChn  = [1:16] + 48;
        if strcmp(date, '20191018')
            fcChn  = [1:15]; 
            pulChn = [1:16] + 16;%line 20
            ppcChn = [1:3,6,10,14:16] + 32;
            vcChn  = [2:5,7,9,11:16] + 48;
        elseif strcmp(date, '20191024') && strcmp(sessionID, '08')
            fcChn  = [1,3:4,6,16]; 
            pulChn = [1:16] + 16;%line
            ppcChn = [1:16] + 32;
            vcChn  = [1:6,11:16] + 48;
        elseif strcmp(date, '20191024') && strcmp(sessionID, '09')
            fcChn  = [1,4:8,11:13,16]; 
            pulChn = [1:16] + 16;%line
            ppcChn = [1:12,14:16] + 32;
            vcChn  = [1:6,10:11,13:16] + 48;
        elseif strcmp(date, '20191025') && strcmp(sessionID, '10')
            fcChn  = [1,4:5,7:8,10,12:13,16]; 
            pulChn = [1,3:12,16] + 16;%line
            ppcChn = [1:12,14:16] + 32;
            vcChn  = [1:6,10,13:16] + 48;
        elseif strcmp(date, '20191025') && strcmp(sessionID, '11')            
            pulChn = [1,3:12,16] + 16;%line
            fcChn  = [1,4:5,7:8,10,12:13,16]; 
            ppcChn = [1:9,14:16] + 32;
            vcChn  = [1:5,13:16] + 48;
        elseif strcmp(date, '20191028') && strcmp(sessionID, '12')            
            pulChn = [1,6,8:13,16] + 16;%line
            fcChn  = [1,3:13,16]; 
            ppcChn = [1:9,14:16] + 32;
            vcChn  = [1:5,13:16] + 48;%weak line
        elseif strcmp(date, '20191029') && strcmp(sessionID, '13')            
            pulChn = [1,3:15,16] + 16;%line 20
            fcChn  = [4:8,10:11,16]; 
            ppcChn = [1:5,9,12,14:16] + 32;
            vcChn  = [1:5,10,13:16] + 48;
        elseif strcmp(date, '20191030') && strcmp(sessionID, '14')            
            pulChn = [1,6:9,12,16] + 16;%line 20
            fcChn  = [1,4:13,16]; 
            ppcChn = [1:5,9,12,14:16] + 32;
            vcChn  = [1:6,10,13:16] + 48; %weak line
        elseif strcmp(date, '20191104') && strcmp(sessionID, '15')            
            pulChn = [1,6:9,12,16] + 16;%line 
            fcChn  = [1,4:10,12:14,16]; 
            ppcChn = [1:4,7,9:10,12,15:16] + 32;
            vcChn  = [1:4,10,13:16] + 48; %weak line
        elseif strcmp(date, '20191104') && strcmp(sessionID, '16')            
            pulChn = [1,3:10,16] + 16;%line 
            fcChn  = [1,5,7:8,12:14,16]; 
            ppcChn = [1:4,7,9:10,12,15:16] + 32;
            vcChn  = [1:9,15:16] + 48; %weak line            
        elseif strcmp(date, '20191105') && strcmp(sessionID, '17')            
            pulChn = [1:10,16] + 16;%line 
            fcChn  = [1,3,5,8,11:13,16]; 
            ppcChn = [1:3,5,9,15:16] + 32;
            vcChn  = [1:4,6,10,13:16] + 48; %weak line
        elseif strcmp(date, '20191106') && strcmp(sessionID, '18')            
            pulChn = [1:16] + 16;%line 
            fcChn  = [1,3:4,9,11:12,14,16]; 
            ppcChn = [1:5,12,15:16] + 32;
            vcChn  = [1:4,10,13:16] + 48; %weak line  
        elseif strcmp(date, '20191106') && strcmp(sessionID, '19')            
            pulChn = [1,3:4,6,10,12:13,15:16] + 16;%line 
            fcChn  = [1,3:4,11:12,14,16]; 
            ppcChn = [1:5,12,15:16] + 32;
            vcChn  = [1:4,10,13:16] + 48; %weak line
        elseif strcmp(date, '20191107') 
            pulChn = [1:12,15:16] + 16;
            fcChn  = [1,5:7,9:10,14:16]; 
            ppcChn = [1:3,5,9,12,15:16] + 32;
            vcChn  = [1:16] + 48;
        elseif strcmp(date, '20191108') && strcmp(sessionID, '21')
            pulChn = [1:11,14:16] + 16; %noisy
            fcChn  = [1:3,6:7,14:16]; %noisy
            ppcChn = [1:5,9:10,12,15:16] + 32; %noisy
            vcChn  = [1:4,6:16] + 48; %noisy
        elseif strcmp(date, '20191108') && strcmp(sessionID, '22')
            pulChn = [1:5,10:11,15:16] + 16;
            fcChn  = [1:3,6:7,15:16];
            ppcChn = [1:5,12,15:16] + 32; 
            vcChn  = [1:4,6:16] + 48;          
        elseif strcmp(date, '20191111') 
            pulChn = [1:5,10:11,15:16] + 16; %noisy
            fcChn  = [1:3,6,16]; %noisy
            ppcChn = [1:5,9:10,12,15:16] + 32; %noisy
            vcChn  = [1:16] + 48; %noisy
        elseif strcmp(date, '20191112') 
            pulChn = [1:5,10:11,15:16] + 16; 
            fcChn  = [1:3,16]; 
            ppcChn = [1:5,12,15:16] + 32; 
            vcChn  = [1:16] + 48;          
        elseif strcmp(date, '20191113') 
            pulChn = [1:4,10:11,15:16] + 16; 
            fcChn  = [1:3,16]; 
            ppcChn = [1:5,7,15:16] + 32; 
            vcChn  = [1:16] + 48;              
        elseif strcmp(date, '20191114') 
            pulChn = [1:2,4,10:11,15:16] + 16; 
            fcChn  = [1:3,16]; 
            ppcChn = [1:5,7:8,15:16] + 32; 
            vcChn  = [1:16] + 48; 
        elseif strcmp(date, '20191115') 
            pulChn = [1:4,10:11,15:16] + 16; 
            fcChn  = [1:2,5:7,12,14,16]; 
            ppcChn = [1:4,15:16] + 32; 
            vcChn  = [1:16] + 48;             
        elseif strcmp(date, '20191118') 
            pulChn = [1:4,10:11,15:16] + 16; 
            fcChn  = [1:2,10,12,16]; 
            ppcChn = [1:2,4,14:16] + 32; 
            vcChn  = [1:16] + 48;              
        elseif strcmp(date, '20191119') 
            pulChn = [1,9:16] + 16; 
            fcChn  = [1:2,10,12,16]; 
            ppcChn = [1:2,15:16] + 32; 
            vcChn  = [1:16] + 48;      
        elseif strcmp(date, '20191120') 
            pulChn = [1:2,4,9:16] + 16; 
            fcChn  = [1:2,5,7,12,16]; 
            ppcChn = [1:2,4:5,9,15:16] + 32; 
            vcChn  = [1:16] + 48;              
        elseif strcmp(date, '20191121') 
            pulChn = [1:4,10:16] + 16; 
            fcChn  = [1:2,8:16]; 
            ppcChn = [1:2,4:5,15:16] + 32; 
            vcChn  = [1:16] + 48;            
        elseif strcmp(date, '20191122') && strcmp(sessionID, '01')
            pulChn = [1:4,10:16] + 16; 
            fcChn  = [1:2,8:16]; 
            ppcChn = [1:2,4:5,9:16] + 32; 
            vcChn  = [1:16] + 48;
        elseif strcmp(date, '20191122') 
            pulChn = [1:4,10:16] + 16; 
            fcChn  = [1:2,8:16]; 
            ppcChn = [1:16] + 32; 
            vcChn  = [1:16] + 48;      
        elseif strcmp(date, '20191125') % new headstages
            pulChn = [1:4,6:10,12:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [1:6,9:16] + 32; 
            vcChn  = [1:16] + 48;  
        elseif strcmp(date, '20191126') 
            pulChn = [1:2,4:6,8,10,12:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [1:6,9,11:16] + 32; 
            vcChn  = [1:16] + 48;                         
        elseif strcmp(date, '20191127') && strcmp(sessionID, '04')
            pulChn = [1:2,4,10:16] + 16; 
            fcChn  = [1:10,12:16]; 
            ppcChn = [1:2,4:6,8,11:16] + 32; 
            vcChn  = [1:16] + 48;             
        elseif strcmp(date, '20191127') && strcmp(sessionID, '01')
            pulChn = [1:2,4,6,10:16] + 16; 
            fcChn  = [1:10,12:16]; 
            ppcChn = [1:3,5:6,11:16] + 32; 
            vcChn  = [1:16] + 48;                  
        elseif strcmp(date, '20191127') && strcmp(sessionType, 'Opto')%opto5-10
            pulChn = [1:6,8,10,12:16] + 16;
            fcChn  = [1:7,9:10,12:16]; 
            ppcChn = [1:6,8:16] + 32;
            vcChn  = [1:16] + 48;
        elseif strcmp(date, '20191128') && strcmp(sessionID, '05')   
            pulChn = [1,4,9:10,12:16] + 16; 
            fcChn  = [1:6,8:9,11:16]; 
            ppcChn = [1:2,5:9,11:16] + 32; 
            vcChn  = [1:16] + 48;
        elseif strcmp(date, '20191128') && strcmp(sessionID, '02') 
            pulChn = [1:2,4,10,12:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [1:2,4:6,8,11:16] + 32; 
            vcChn  = [1:8,10:16] + 48;
        elseif strcmp(date, '20191129') && strcmp(sessionID, '06')
            pulChn = [1,4,10:12,15:16] + 16; 
            fcChn  = [1:9,11:16]; 
            ppcChn = [1:3,5,11:16] + 32; 
            vcChn  = [1,3:14,16] + 48;
        elseif strcmp(date, '20191129')    
            pulChn = [1:7,9,12:16] + 16; 
            fcChn  = [1:6,8:9,11:16]; 
            ppcChn = [1:3,5,7,10:16] + 32; 
            vcChn  = [1:14,16] + 48;
        elseif strcmp(date, '20191129') && strcmp(sessionType, 'Opto')%opto11-16
            pulChn = [1,3:5,9:16] + 16;
            fcChn  = [1:9,11:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;
        elseif strcmp(date, '20191203') 
            pulChn = [1,4,9:12,15:16] + 16; 
            fcChn  = [1:7,9:11,13:16]; 
            ppcChn = [1:3,10:16] + 32; 
            vcChn  = [1,3:6,8,10:12,14,16] + 48; 
        elseif strcmp(date, '20191204') 
            pulChn = [1,4,10:12,15:16] + 16; 
            fcChn  = [1:7,9:16]; 
            ppcChn = [1:3,5,9:16] + 32; 
            vcChn  = [1:6,8:14,16] + 48;     
        elseif strcmp(date, '20191205') 
            pulChn = [1:16] + 16; 
            fcChn  = [1:7,9:16]; 
            ppcChn = [1:3,5,9:16] + 32; 
            vcChn  = [1,3:6,9:14,16] + 48; 
        elseif strcmp(date, '20191206') 
            pulChn = [1:16] + 16;
            fcChn  = [1,3:7,10:11,13:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1,3:6,9,11:13,16] + 48;
        elseif strcmp(date, '20191209') 
            pulChn = [1:16] + 16;
            fcChn  = [1,3:7,10:11,13,16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:5,7:16] + 48;   
        elseif strcmp(date, '20191211') 
            pulChn = [1:16] + 16;
            fcChn  = [1,4,6:7,10,13,16]; 
            ppcChn = [1,3:16] + 32;
            vcChn  = [1:5,7:8,10:16] + 48;  
        elseif strcmp(date, '20191212') 
            pulChn = [1:16] + 16;
            fcChn  = [1,4:7,9,11,13,16]; 
            ppcChn = [1,3:16] + 32;
            vcChn  = [1:5,7:8,10:16] + 48;
        elseif strcmp(date, '20191213') % both sessions
            pulChn = [1:10,12,14:16] + 16;
            fcChn  = [1,4:7,10,16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:5,9,11:16] + 48;
        elseif strcmp(date, '20191216') && strcmp(sessionID, '11')
            pulChn = [1:10,12,14:16] + 16;
            fcChn  = [1,4:8,10:11,16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:2,4:7,9,11:16] + 48;    
        elseif strcmp(date, '20191216') && strcmp(sessionID, '09')
            pulChn = [1:16] + 16;
            fcChn  = [1:14,16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1,5,7,9,11:16] + 48; 
        elseif strcmp(date, '20191217')
            pulChn = [1:9,12,14:16] + 16;
            fcChn  = [1:14,16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1,5,9,11:13,16] + 48;             
        elseif strcmp(date, '20191218')
            pulChn = [1:9,12:16] + 16;
            fcChn  = [1:14,16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:8,10:16] + 48;
        elseif strcmp(date, '20191219')
            pulChn = [1:9,12:16] + 16;
            fcChn  = [1:13,16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;    % noisy spikes
        elseif strcmp(date, '20191221')
            pulChn = [1:9,12:16] + 16;
            fcChn  = [1:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48; 
        elseif strcmp(date, '20191222')
            pulChn = [1:9,12:16] + 16;
            fcChn  = [1:14,16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48; 
        elseif strcmp(date, '20191222') && strcmp(sessionID, '15')
            pulChn = [1:9,12:16] + 16;
            fcChn  = [1:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:11,13:15] + 48; 
        elseif strcmp(date, '20191223')
            pulChn = [1:10,12:16] + 16;
            fcChn  = [1:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:2,4:10,12:15] + 48;             
        elseif strcmp(date, '20191228')
            pulChn = [1:16] + 16; %noisy
            fcChn  = [1:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;       
        elseif strcmp(date, '20191229')
            pulChn = [1:16] + 16; %noisy
            fcChn  = [1:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;     
        elseif strcmp(date, '20200102')
            pulChn = [1:10,12:16] + 16; %noisy
            fcChn  = [1:9,11:13,15:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;             
        elseif strcmp(date, '20200103') && strcmp(sessionID, '17')
            pulChn = [1:16] + 16; 
            fcChn  = [1:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:6,8:16] + 48;
        elseif strcmp(date, '20200103')
            pulChn = [1:10,12:16] + 16; %noisy
            fcChn  = [1:10,12:13,15:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:6,8:16] + 48; % -- above all filled 1/4/2020 
        elseif strcmp(date, '20200109') && strcmp(sessionID, '18')
            pulChn = [1:7,9:10,12:16] + 16;
            fcChn  = [1:12,14:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:6,8:9,11:16] + 48;    
        elseif strcmp(date, '20200108') || strcmp(date, '20200109')
            pulChn = [1:10,12:16] + 16;
            fcChn  = [1:12,14:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:6,8:16] + 48;            
        elseif strcmp(date, '20200110')
            pulChn = [1:7,9:10,12:16] + 16;
            fcChn  = [1:5,9:12,14:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;  
        elseif strcmp(date, '20200113') && strcmp(sessionID, '13')
            pulChn = [1:16] + 16;
            fcChn  = [1:6,9:12,14:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48; 
        elseif strcmp(date, '20200113')
            pulChn = [1:7,9:10,12:16] + 16;
            fcChn  = [1:6,9:12,14:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48; 
        elseif strcmp(date, '20200114')
            pulChn = [1:7,9:10,12:16] + 16;
            fcChn  = [1:5,9:12,14:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:8,10:16] + 48;
        elseif strcmp(date, '20200115') && strcmp(sessionID, '08')
            pulChn = [1:8,10,12:16] + 16;
            fcChn  = [1:5,9:12,14:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:8,10:16] + 48;
        elseif strcmp(date, '20200115') 
            pulChn = [1:5,9:10,12:16] + 16;
            fcChn  = [1:4,14:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:8,10:16] + 48;    
        elseif strcmp(date, '20200116') 
            pulChn = [1:6,10,12:16] + 16;
            fcChn  = [1:5,14:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:2,4:8,10:16] + 48;     
        elseif strcmp(date, '20200117') && strcmp(sessionID, '10')
            pulChn = [1:10,13:16] + 16;
            fcChn  = [1:3,5,16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:2,4:8,10:16] + 48; 
        elseif strcmp(date, '20200117')
            pulChn = [1:10,12:16] + 16;
            fcChn  = [1:3,5,16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;    
        elseif strcmp(date, '20200121') && strcmp(sessionID, '11')
            pulChn = [1:10,13:16] + 16;
            fcChn  = [1:3,16]; 
            ppcChn = [1:6,8:14,16] + 32;
            vcChn  = [1:14,16] + 48;     
        elseif strcmp(date, '20200121')
            pulChn = [1:10,13:16] + 16;
            fcChn  = [1:3,16]; 
            ppcChn = [1:7,9:16] + 32;
            vcChn  = [1:14,16] + 48;   
        elseif strcmp(date, '20200122') && strcmp(sessionID, '12')
            pulChn = [1:10,12:16] + 16;
            fcChn  = [1:5,16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;     
        elseif strcmp(date, '20200122')
            pulChn = [1:10,12:16] + 16;
            fcChn  = [1:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;      
        elseif strcmp(date, '20200123')
            pulChn = [1:10,12:16] + 16;
            fcChn  = [1:14,16]; 
            ppcChn = [1:7,9:16] + 32;
            vcChn  = [1:13,15:16] + 48;
        elseif strcmp(date, '20200124') %both
            pulChn = [1:10,12:16] + 16;
            fcChn  = [1:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;    
        elseif strcmp(date, '20200127') && strcmp(sessionID, '21')
            pulChn = [1:10,12:16] + 16;
            fcChn  = [1:16]; 
            ppcChn = [1:3,6:16] + 32;
            vcChn  = [1:2,4:16] + 48;   
        elseif strcmp(date, '20200127')
            pulChn = [1:10,12:16] + 16;
            fcChn  = [1:16]; 
            ppcChn = [1:5,6:13,15:16] + 32;
            vcChn  = [1:2,4:14,16] + 48;   
        elseif strcmp(date, '20200128')
            pulChn = [1:10,12:16] + 16;
            fcChn  = [1:16]; 
            ppcChn = [1:3,6:7,9:16] + 32;
            vcChn  = [1:2,4:11,13,16] + 48;   
        elseif strcmp(date, '20200129') && strcmp(sessionID, '22')
            pulChn = [1:10,12:16] + 16;
            fcChn  = [1:15]; 
            ppcChn = [1:5,7:10,12,14:16] + 32;
            vcChn  = [1:11,13:14,16] + 48; 
        elseif strcmp(date, '20200129')
            pulChn = [1:10,12:16] + 16;
            fcChn  = [1:15]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:11,13:14,16] + 48;     
        elseif strcmp(date, '20200131')
            pulChn = [1:10,12:16] + 16;
            fcChn  = [1:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:2,4:16] + 48;   % -- above all filled on 2/1/2020
        elseif strcmp(date, '20200131')
            pulChn = [1:10,12:16] + 16;
            fcChn  = [1:16]; 
            ppcChn = [1:16] + 32;
            vcChn  = [1:2,4:16] + 48;   
        elseif strcmp(date, '20200203')
            pulChn = [1:10,12:16] + 16;
            fcChn  = [1:16]; 
            ppcChn = [1:10,12:16] + 32;
            vcChn  = [1:2,4:16] + 48;  
        elseif strcmp(date, '20200204') % PFC noisy
            pulChn = [1:10,12:16] + 16;
            fcChn  = [1:16]; 
            ppcChn = [1:5,7:12,14:16] + 32;
            vcChn  = [1:2,4:16] + 48;  
        elseif strcmp(date, '20200205') % PFC noisy
            pulChn = [1:10,12:16] + 16;
            fcChn  = [1:16]; 
            ppcChn = [1:5,7:12,14:16] + 32;
            vcChn  = [1:2,4:16] + 48;  
        elseif strcmp(date, '20200206')
            pulChn = [1:10,12:16] + 16;
            fcChn  = [1:16]; 
            ppcChn = [1:4,7,9:10,12:16] + 32;
            vcChn  = [1:16] + 48; 
        elseif strcmp(date, '20200207')
            pulChn = [1:9,12:16] + 16;
            fcChn  = [1:16]; 
            ppcChn = [1:7,9:10,12:16] + 32;
            vcChn  = [1:16] + 48;    
        elseif strcmp(date, '20200213') %level8 both
            pulChn = [1:10,12:16] + 16;
            fcChn  = [1:7,9:10,12,14:16]; 
            ppcChn = [1:4,7,9,12:13,15:16] + 32;
            vcChn  = [1:2,4:6,8:16] + 48;              
        elseif strcmp(date, '20200214') 
            pulChn = [1:5,7:9,12:16] + 16;
            fcChn  = [1:7,9,12,15:16]; 
            ppcChn = [1:5,7,9,12:16] + 32;
            vcChn  = [1:2,4:6,8:16] + 48; % -- above all filled 2/17/2020
        elseif strcmp(date, '20200217')
            pulChn = [1:4,7:10,12:16] + 16; 
            fcChn  = [1:6,9,12:13,15:16]; 
            ppcChn = [1:4,7,9,12:16] + 32;
            vcChn  = [1:4,6,8:16] + 48;     
        elseif strcmp(date, '20200218')
            pulChn = [1:4,7:10,12:16] + 16; 
            fcChn  = [1:6,12:13,16]; 
            ppcChn = [1:4,7,9,12:16] + 32;
            vcChn  = [1:4,8:16] + 48;
        elseif strcmp(date, '20200219')
            pulChn = [1,3:5,7:10,12:16] + 16; 
            fcChn  = [1:5,12:13,16]; 
            ppcChn = [1:4,7,9,12:16] + 32;
            vcChn  = [1:4,8:12,15:16] + 48;
        elseif strcmp(date, '20200220') && strcmp(sessionID, '01')% both level9
            pulChn = [1,3:10,12:16] + 16; 
            fcChn  = [1,3,5:6,9,12:13,16]; 
            ppcChn = [1:4,7,12:14,16] + 32;
            vcChn  = [1:4,9,11:12,15:16] + 48; 
        elseif strcmp(date, '20200220') 
            pulChn = [1,3:10,12:14,16] + 16; 
            fcChn  = [1,3,5:6,9,13,16]; 
            ppcChn = [1:4,7,12:13,15:16] + 32;
            vcChn  = [1:4,8:9,11:12,16] + 48;
        elseif strcmp(date, '20200224') && strcmp(sessionID, '03')
            pulChn = [1,3:10,12:14,16] + 16; 
            fcChn  = [1,3,5:6,9,13,16]; 
            ppcChn = [1:7,12:13,16] + 32;
            vcChn  = [1:5,8,11:12,16] + 48;
        elseif strcmp(date, '20200224') 
            pulChn = [1,3:10,12:14,16] + 16; 
            fcChn  = [1,3,9,13,16]; 
            ppcChn = [1:4,7,12:16] + 32;
            vcChn  = [1:4,11:13,16] + 48;
        elseif strcmp(date, '20200225') 
            pulChn = [1,3:9,12:16] + 16; 
            fcChn  = [1,3,9,13,16]; 
            ppcChn = [1:4,7,9,12:16] + 32;
            vcChn  = [1:4,11:13,16] + 48;
        elseif strcmp(date, '20200302') 
            pulChn = [1,3:10,12:14,16] + 16; 
            fcChn  = [1,16]; 
            ppcChn = [1:5,7,12:13,15:16] + 32;
            vcChn  = [1:2,4,11:12,14,16] + 48;
        elseif strcmp(date, '20200303') 
            pulChn = [1,3:10,12:14,16] + 16; 
            fcChn  = [1,16]; 
            ppcChn = [1:4,7,12:16] + 32;
            vcChn  = [1:2,4,11:12,16] + 48;
        elseif strcmp(date, '20200304') && strcmp(sessionID, '08')
            pulChn = [1,4:5,7:14,16] + 16; 
            fcChn  = [1:12,14:16]; 
            ppcChn = [1:4,7,12:13,16] + 32;
            vcChn  = [1:2,8,11:13,16] + 48;
        elseif strcmp(date, '20200304') 
            pulChn = [1,4,7:10,12:14,16] + 16; 
            fcChn  = [1:12,14:16]; 
            ppcChn = [1:4,12,14:16] + 32;
            vcChn  = [1:2,4,11:12,16] + 48;
        elseif strcmp(date, '20200305') 
            pulChn = [1,4,7:10,12:14,16] + 16; 
            fcChn  = [1:9,11:12,14:16]; 
            ppcChn = [1:4,12:14,16] + 32;
            vcChn  = [1:4,9,11:13,16] + 48;    
        elseif strcmp(date, '20200306') 
            pulChn = [1,8:10,12:14,16] + 16; 
            fcChn  = [1:9,12,14:16]; 
            ppcChn = [1:4,7,12,14,16] + 32;
            vcChn  = [1:4,6:13,16] + 48; % -- above all filled 3/6/2020
        end
        %validChn = {fcChn; pulChn; ppcChn; vcChn}; % end of 0181
        
        % Check if pul channel is annatomically valid (added 3/26/2021)
        validSUMask = getValidSUMask(animalCode,'LPl',pulChn-16);
        newpulChn = pulChn(validSUMask);
        validChn = {fcChn; newpulChn; ppcChn; vcChn};
%%            
    case '0171'        
        fcChn  = [1:16]; %'20190304'
        pulChn = [1:12,14:16] + 16;
        ppcChn = [1:16] + 32;
        vcChn  = [1:16] + 48;
        if strcmp(date, '20190305')
            fcChn  = [1,3:10,12:16];
            ppcChn = [2:7,9:16] + 32;
            vcChn  = [2,5:8,11:13,15:16] + 48;
        elseif strcmp(date, '20190306') || strcmp(date, '20190307') || strcmp(date, '20190308') || strcmp(date, '20190312')  
            %all use default setting
        elseif strcmp(date, '20190313') || strcmp(date, '20190314') || strcmp(date, '20190315')
            ppcChn = [1:12,14:16] + 32;
        elseif strcmp(date, '20190318') || strcmp(date, '20190319')
            ppcChn = [1:10,12,14:16] + 32;
            vcChn  = [1:4,6:16] + 48;
        elseif strcmp(date, '20190320')
            pulChn = [1:7,9:16] + 16;
            ppcChn = [1:10,12,14:16] + 32;
        elseif strcmp(date, '20190321')
            pulChn = [1:16] + 16;
            ppcChn = [1:4,6:10,12,14:16] + 32;            
        elseif strcmp(date, '20190327')
            pulChn = [1:7,10:16] + 16;
            ppcChn = [1,3:4,6:10,12,14,16] + 32;
        elseif strcmp(date, '20190328') || strcmp(date, '20190329') || strcmp(date, '20190401') || strcmp(date, '20190402')
            pulChn = [1:6,10:16] + 16; %PPC is good again
            ppcChn = [1:2,4:16] + 32;
        elseif strcmp(date, '20190403')
            pulChn = [1:7,10:16] + 16;
            ppcChn = [1:2,4:8,10:16] + 32;
        elseif strcmp(date, '20190405')
            pulChn = [1:7,10:16] + 16;
            vcChn  = [1:2,5:16] + 48;
        elseif strcmp(date, '20190408')
            pulChn = [1:7,10:16] + 16;
        elseif strcmp(date, '20190409')
            pulChn = [1:6,10:16] + 16;
            ppcChn = [1:14,16] + 32;
            vcChn  = [1:2,4:16] + 48;
        elseif strcmp(date, '20190419')
            pulChn = [1:7,10,13:15] + 16;
            ppcChn = [1:9,13:16] + 32;
            vcChn  = [1:2,4:16] + 48;
        elseif strcmp(date, '20190422') || strcmp(date, '20190423') || strcmp(date, '20190424') || strcmp(date, '20190426')
            fcChn  = [1:12,14:16]; 
            pulChn = [1:16] + 16;
            ppcChn = [1:9,13:16] + 32;
            vcChn  = [1:2,4:16] + 48;
        elseif strcmp(date, '20190425') % 0426 is above
            fcChn  = [1:12,14:16]; 
            pulChn = [1:16] + 16;
            ppcChn = [1:8,10:14] + 32;
            vcChn  = [1:2,4:16] + 48;
        elseif strcmp(date, '20190427') || strcmp(date, '20190429') || strcmp(date, '20190430')
            fcChn  = [1:12,14:16]; 
            pulChn = [1:7,10:16] + 16;
            ppcChn = [1:10,12:16] + 32;
            vcChn  = [1:2,4:7,9:16] + 48;    
        elseif strcmp(date, '20190503')
            fcChn  = [1:12,14:16]; 
            pulChn = [1:7,10:16] + 16;
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;   
        elseif strcmp(date, '20190507')
            fcChn  = [1:9,11:12,14:16]; 
            pulChn = [1:7,10:16] + 16;
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;  
        elseif strcmp(date, '20190510')
            fcChn  = [1:9,11:12,14:16]; 
            pulChn = [1:7,10:16] + 16;
            ppcChn = [1:7,14:16] + 32;
            vcChn  = [1:7,9:16] + 48;  
        elseif strcmp(date, '20190513')
            fcChn  = [1:12,14:16]; 
            pulChn = [1:7,10:16] + 16;
            ppcChn = [1:10,12:16] + 32;
            vcChn  = [1:7,9:16] + 48;
        elseif strcmp(date, '20190515')
            fcChn  = [1:12,14:16]; 
            pulChn = [1:3,6:7,10:16] + 16;
            ppcChn = [1:10,12:16] + 32;
            vcChn  = [1:7,10:16] + 48;
        elseif strcmp(date, '20190516')
            fcChn  = [1:12,14:16]; 
            pulChn = [1:7,10:16] + 16;
            ppcChn = [1:10,13:16] + 32;
            vcChn  = [1:4,6:7,9,11:16] + 48;     
        elseif strcmp(date, '20190517')
            fcChn  = [1:8,10:16]; 
            pulChn = [1:7,10:16] + 16;
            ppcChn = [1:2,6:8,10:11,14:16] + 32;
            vcChn  = [1:4,6:7,9,11:16] + 48;   
        elseif strcmp(date, '20190521')
            fcChn  = [1:7,9:12,14:16]; 
            pulChn = [1:7,10:16] + 16;
            ppcChn = [1:8,13:14,16] + 32;
            vcChn  = [1:2,4,6:7,11:16] + 48;          
        elseif strcmp(date, '20190523')
            fcChn  = [1:7,9:12,14:16]; 
            pulChn = [1:6,10:11,14:16] + 16;
            ppcChn = [1:8,13:16] + 32;
            vcChn  = [1:4,6:7,11:16] + 48;
        elseif strcmp(date, '20190524')
            fcChn  = [1:7,9:12,14:16]; 
            pulChn = [1:7,10:16] + 16;
            ppcChn = [1:8,13:16] + 32;
            vcChn  = [1:4,6:7,11:16] + 48;            
        elseif strcmp(date, '20190527')
            fcChn  = [1:7,9:12,14:16]; 
            pulChn = [1:6,10:11,14:16] + 16;
            ppcChn = [1:2,4:8,13,15:16] + 32;
            vcChn  = [1:4,6:7,11:12,14:16] + 48;             
        elseif strcmp(date, '20190528')
            fcChn  = [1:2,5:16]; 
            pulChn = [1:4,6,10:12,14:16] + 16;
            ppcChn = [1:4,6:8,13:15] + 32;
            vcChn  = [1:4,6:7,11:16] + 48;
        elseif strcmp(date, '20190529')
            fcChn  = [1:7,9:12,14:16]; 
            pulChn = [1:4,6,10:16] + 16;
            ppcChn = [1:14,16] + 32;
            vcChn  = [1:4,6:7,11:16] + 48;
        elseif strcmp(date, '20190530')
            fcChn  = [1:7,9:12,14:16]; 
            pulChn = [1:4,6,10:16] + 16;
            ppcChn = [1:4,6:10,13,16] + 32;
            vcChn  = [1:4,6,11:12,14:16] + 48;             
        elseif strcmp(date, '20190531')
            fcChn  = [1:7,9:12,14:16]; 
            pulChn = [1:4,6,10:11,14:16] + 16;
            ppcChn = [1:3,5,7,9:11,13,15:16] + 32;
            vcChn  = [1:6,10:12,14:16] + 48;            
        elseif strcmp(date, '20190603')
            fcChn  = [1:7,9:12,14:16]; 
            pulChn = [1:4,6,10:11,14:16] + 16;
            ppcChn = [1:16] + 32;
            vcChn  = [1:4,12,14:16] + 48;                 
        elseif strcmp(date, '20190604')
            fcChn  = [1:7,9:12,14:16]; 
            pulChn = [1:4,6,14:16] + 16;
            ppcChn = [1:9,11:16] + 32;
            vcChn  = [1:4,12,14:16] + 48;
        elseif strcmp(date, '20190607')
            fcChn  = [1:7,9:12,14:16]; 
            pulChn = [1:4,6,11,14:16] + 16;
            ppcChn = [1:9,11:16] + 32;
            vcChn  = [1:4,12,14:16] + 48;
        end
        %validChn = {fcChn; pulChn; ppcChn; vcChn}; % end of 0171
    
        % Check if pul channel is annatomically valid (added 3/26/2021)
        validSUMask = getValidSUMask(animalCode,'LPl',pulChn-16);
        newpulChn = pulChn(validSUMask);
        validChn = {fcChn; newpulChn; ppcChn; vcChn};
    %%
    case '0173'
        fcChn  = [1:13,15:16];
        pulChn = [1:9,13:16] + 16; %excluding (17,21,22)+
        ppcChn = [1:16] + 32;
        vcChn  = [1:16] + 48;
        if strcmp(date, '20190104')
            fcChn  = [1:12,14:16];
            pulChn = [1:3,5:10,13:16] + 16;
        elseif strcmp(date, '20190105') %
            vcChn  = [1:12,14:16] + 48;
        elseif strcmp(date, '20190107') %
            fcChn  = [1:16];
            vcChn  = [1:9,11:12,14,16] + 48;
        elseif strcmp(date, '20190108') %
            fcChn  = [1:10,12:15];
            vcChn  = [1:9,11:12,14,16] + 48;
        elseif strcmp(date, '20190109')
            fcChn  = [1:3,5,7:16];
            pulChn = [1:3,5:16] + 16;
            ppcChn = [1:3,5:16] + 32;
            vcChn  = [1:9,11:12,14,16] + 48;
        elseif strcmp(date, '20190110')
            fcChn  = [1:2,4:16];
            pulChn = [1:15] + 16;
            ppcChn = [1:2,5:16] + 32;
            vcChn  = [1:5,9,11:12,14,16] + 48;
        elseif strcmp(date, '20190111')
            fcChn  = [1:2,4:16];
            pulChn = [1:6,8:11,13:14] + 16;
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;
        elseif strcmp(sessionID, '10')
            fcChn  = [1:5,7:16];
            pulChn = [1:6,8:11,13:14] + 16;
            ppcChn = [1:16] + 32;
            vcChn  = [1:14,16] + 48;
        elseif strcmp(sessionID, '11')
            fcChn  = [1:5,7:16];
            pulChn = [1,3:5,8] + 16;
            ppcChn = [1:14,16] + 32;
            vcChn  = [1:2,4:14,16] + 48;
        elseif strcmp(date, '20190113')
            fcChn  = [1:5,7:16];
            pulChn = [1,3,5,13] + 16;
            ppcChn = [1:14,16] + 32;
            vcChn  = [1:2,4:14,16] + 48;
        elseif strcmp(date, '20190114') %movement artifact
            fcChn  = [1:5,7:16];
            pulChn = [1,5,8] + 16;
            ppcChn = [1:14,16] + 32;
            vcChn  = [1,4,6:10,13:14,16] + 48;
        elseif strcmp(date, '20190115') %movement artifact
            fcChn  = [1:5,8:16];
            pulChn = [1,5,8] + 16;
            ppcChn = [1:10] + 32;
            vcChn  = [1,4,9:10,13:14,16] + 48;
        elseif strcmp(date, '20190116') %movement artifact
            fcChn  = [1:5,8:16]; 
            pulChn = [1,5,8] + 16;
            ppcChn = [1:10] + 32;
            vcChn  = [1,4,16] + 48;    
        elseif strcmp(date, '20190117') %movement artifact
            fcChn  = [1:5,8:16]; 
            pulChn = [1,5] + 16;
            ppcChn = [1:16] + 32;
            vcChn  = [1,4,16] + 48;     
        elseif strcmp(date, '20190118')
            fcChn  = [1:5,8:16]; 
            pulChn = [5:8,11,13:15] + 16;
        elseif strcmp(date, '20190121') 
            pulChn = [4:8,11,13:15] + 16;
        elseif strcmp(date, '20190122') 
            pulChn = [4:8,11,13:15] + 16;        
            ppcChn = [1:9,11:16] + 32;
            vcChn  = [1:9,11,13:16] + 48;
        elseif strcmp(date, '20190123') % signal looks saturated
            pulChn = [4:8,11,13:15] + 16;           
            vcChn  = [1:9,11,13:16] + 48;
        elseif strcmp(date, '20190124') % signal looks saturated
            fcChn  = [1:5,8:16]; 
            pulChn = [4:8,11,13:15] + 16;           
            vcChn  = [1:9,11,13:16] + 48;
        elseif strcmp(date, '20190126') % signal looks saturated
            fcChn  = [1:9,15:16]; 
            pulChn = [4:8,13:15] + 16;
            ppcChn = [1:8,11:16] + 32;
            vcChn  = [1:4,6:7] + 48;            
        elseif strcmp(date, '20190128')
            fcChn  = [1:6,8:16];
            pulChn = [2:8,12:16] + 16;
        elseif strcmp(date, '20190129')
            fcChn  = [1:6,8:16];
            pulChn = [2:8,12:16] + 16;
            ppcChn = [1:13,15:16] + 32;
        elseif strcmp(date, '20190130') || strcmp(date, '20190131')
            fcChn  = [1:16];
            pulChn = [2:8,12:16] + 16;
        elseif strcmp(date, '20190201')
            fcChn  = [1:16];
            pulChn = [2:8,13:16] + 16;
            ppcChn = [1:9,11:16] + 32;
        elseif strcmp(date, '20190204')
            fcChn  = [1:16];
            pulChn = [2:8,13:16] + 16;
            ppcChn = [1:10] + 32;
            vcChn  = [1:8] + 48;
        elseif strcmp(date, '20190205')
            fcChn  = [1:15];
            pulChn = [2:8,13:16] + 16;
            ppcChn = [1:10] + 32;
            vcChn  = [1:8] + 48;
        elseif strcmp(date, '20190206')
            fcChn  = [1:4,6:16];
            pulChn = [2:15] + 16;
            ppcChn = [1:10] + 32;
            vcChn  = [1:8] + 48; 
        elseif strcmp(date, '20190207')
            fcChn  = [1:4,6:14];
            pulChn = [2:11,13:15] + 16;
            ppcChn = [1:14] + 32;
            vcChn  = [1:8] + 48;  
        elseif strcmp(date, '20190208')
            fcChn  = [1:14];
            pulChn = [2:8,13:15] + 16;
            ppcChn = [1:8,13:14] + 32;
            vcChn  = [1:9] + 48; 
        elseif strcmp(date, '20190212')
            fcChn  = [2:4,6:16];
            pulChn = [2:7,10:13,15] + 16;
            ppcChn = [1:10] + 32;
            vcChn  = [1:8] + 48;
        elseif strcmp(date, '20190213')
            fcChn  = [1:12,14:16];
            pulChn = [3:8,11:16] + 16;
            vcChn  = [1:8] + 48;
        elseif strcmp(date, '20190215')
            fcChn  = [1:12,14:16];
            pulChn = [2:8,10,12:16] + 16;
            ppcChn = [1:16] + 32;
            vcChn  = [1:8] + 48;
        elseif strcmp(date, '20190218') || strcmp(date, '20190219') || strcmp(date, '20190220')
            fcChn  = [1:12,14:16];
            pulChn = [2:8,10:16] + 16;
            ppcChn = [1:16] + 32;
            vcChn  = [1:8] + 48;            
        elseif strcmp(date, '20190221')
            fcChn  = [1:16];
            pulChn = [2:8,10:16] + 16;
            ppcChn = [1:9,11:15] + 32;
            vcChn  = [1:4,7:8] + 48;
        elseif strcmp(date, '20190222')
            pulChn = [1:5,7:8,13:15] + 16;
            ppcChn = [1:8,11:15] + 32;
            vcChn  = [1:4,7:8] + 48; 
        elseif strcmp(date, '20190227') || strcmp(date, '20190228')
            fcChn  = [1:10,12:16];
            pulChn = [2:8,14:15] + 16;
            ppcChn = [3:8,11:14] + 32;
            vcChn  = [1:8] + 48; 
        elseif strcmp(date, '20190301')
            fcChn  = [1:10,12:16];
            pulChn = [2:8,10:15] + 16;
            ppcChn = [1:9] + 32;
            vcChn  = [1:8] + 48;
        elseif strcmp(date, '20190304')
            fcChn  = [2:16];
            pulChn = [4:8,11:15] + 16;
            ppcChn = [1,3:16] + 32;
            vcChn  = [1:8,16] + 48;
        elseif strcmp(date, '20190305')
            fcChn  = [2:16];
            pulChn = [4:8,11:15] + 16;
            ppcChn = [1,3:16] + 32;
            vcChn  = [1:8,16] + 48;        
        elseif strcmp(date, '20190307')
            fcChn  = [1:16];
            pulChn = [2:8,10:15] + 16;
            ppcChn = [1:7,9:16] + 32;
            vcChn  = [1:8,10:16] + 48; 
        elseif strcmp(date, '20190312')
            fcChn  = [1:16];
            pulChn = [2:7,13:15] + 16;
            ppcChn = [1:12,14:16] + 32;
            vcChn  = [1:16] + 48;  
        elseif strcmp(date, '20190314')
            fcChn  = [2:6,8:16];
            pulChn = [2:8,10:16] + 16;
            ppcChn = [3:16] + 32;
            vcChn  = [1:11] + 48;  
        elseif strcmp(date, '20190315')||strcmp(date, '20190318')||strcmp(date, '20190319')
            fcChn  = [1:16];
            pulChn = [2:7,11:13,15] + 16;
            ppcChn = [1:10,12,14:15] + 32;
            vcChn  = [1:16] + 48; 
        elseif strcmp(date, '20190320')
            fcChn  = [1:16];
            pulChn = [2:7,10:16] + 16;
            ppcChn = [1:4,6:8,10,12,14:16] + 32;
            vcChn  = [1,3:16] + 48; 
        elseif strcmp(date, '20190326')
            fcChn  = [1:16];
            pulChn = [2:7,10:16] + 16;
            ppcChn = [1:3,6:8,10,12,14,16] + 32;
            vcChn  = [1:16] + 48; 
        elseif strcmp(date, '20190327')
            fcChn  = [1:10,12:16];
            pulChn = [2:7,10:16] + 16;
            ppcChn = [1:3,6:8,10,12,14,16] + 32;
            vcChn  = [1:16] + 48;
        elseif strcmp(date, '20190328') || strcmp(date, '20190329')
            fcChn  = [1:10,12:16];
            pulChn = [2:7,10:12,14:16] + 16;
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;
        elseif strcmp(date, '20190401')
            fcChn  = [1:16];
            pulChn = [2:5,10:16] + 16;
            ppcChn = [1:16] + 32;
            vcChn  = [1:16] + 48;
        elseif strcmp(date, '20190402') || strcmp(date, '20190403')
            fcChn  = [1:16];
            pulChn = [2:7,10:12,14:16] + 16;
            ppcChn = [1:16] + 32;
            vcChn  = [1:10,16] + 48;    
        elseif strcmp(date, '20190405')
            fcChn  = [1:14,16];
            pulChn = [2:7,11:15] + 16;
            ppcChn = [1:16] + 32;
            vcChn  = [1:10,16] + 48; 
        elseif strcmp(date, '20190408')
            fcChn  = [1:16];
            pulChn = [2:7,11:15] + 16;
            ppcChn = [1:16] + 32;
            vcChn  = [1:10,16] + 48; 
        elseif strcmp(date, '20190409') || strcmp(date, '20190410') || strcmp(date, '20190411')
            fcChn  = [1:7,9:16];
            pulChn = [2:7,11:15] + 16;
            ppcChn = [1:16] + 32;
            vcChn  = [1:3,5:10,16] + 48;
        elseif strcmp(date, '20190419')
            fcChn  = [1:7,9:10,12:16];
            pulChn = [2:7,11,13:15] + 16;
            ppcChn = [1:9,15:16] + 32;
            vcChn  = [1:10] + 48;            
        elseif strcmp(date, '20190423')
            fcChn  = [1:7,9:16];
            pulChn = [2:7,11,13:15] + 16;
            ppcChn = [1:9,15:16] + 32;
            vcChn  = [1:10,16] + 48;                      
        elseif strcmp(date, '20190424')
            fcChn  = [1:7,9:16];
            pulChn = [2:7,13:15] + 16;
            ppcChn = [1:9,15:16] + 32;
            vcChn  = [1:9,16] + 48;  
        elseif strcmp(date, '20190425')
            fcChn  = [1:7,9:16];
            pulChn = [2:8,11:15] + 16;
            ppcChn = [1:9,15:16] + 32;
            vcChn  = [1:9] + 48;              
        elseif strcmp(date, '20190426') || strcmp(date, '20190429')
            fcChn  = [1:7,9:16];
            pulChn = [2:8,10:15] + 16;
            ppcChn = [1:16] + 32;
            vcChn  = [1:9] + 48;
        elseif strcmp(date, '20190503') || strcmp(date, '20190507')
            fcChn  = [1:7,9:14,16];
            pulChn = [1:8,10:15] + 16;
            ppcChn = [1:13,15:16] + 32;
            vcChn  = [1:8,10] + 48;
        elseif strcmp(date, '20190510')
            fcChn  = [1:7,9:14,16];
            pulChn = [1:2,5:7,9:15] + 16;
            ppcChn = [1:13,15:16] + 32;
            vcChn  = [1:4,6:10] + 48;            
        elseif strcmp(date, '20190513')
            fcChn  = [1:7,9:14,16];
            pulChn = [2,10:12,14:16] + 16;
            ppcChn = [1:9,12:13,15:16] + 32;
            vcChn  = [1:6,8,16] + 48;
        elseif strcmp(date, '20190516')
            fcChn  = [1:7,9:14,16];
            pulChn = [1:7,9:16] + 16;
            ppcChn = [1:9,12:13,15:16] + 32;
            vcChn  = [1:8,10,15:16] + 48;              
        elseif strcmp(date, '20190517')
            fcChn  = [1:7,9:14,16];
            pulChn = [1:7,10:16] + 16;
            ppcChn = [1:11,15:16] + 32;
            vcChn  = [1:8,10,15:16] + 48;
        elseif strcmp(date, '20190521')
            fcChn  = [1:7,9:14,16];
            pulChn = [2,4:7,13:15] + 16;
            ppcChn = [1:9,11:13,15] + 32;
            vcChn  = [1:6,8,16] + 48;   
        elseif strcmp(date, '20190523')
            fcChn  = [1:7,9:14,16];
            pulChn = [1:7,11:16] + 16;
            ppcChn = [1:9,11:13,15] + 32;
            vcChn  = [1:7,15:16] + 48;  
        elseif strcmp(date, '20190527')
            fcChn  = [1:7,9:14,16];
            pulChn = [1:7,11:16] + 16;
            ppcChn = [1:7,12:13,15] + 32;
            vcChn  = [1:4,6:7,15:16] + 48;        
        elseif strcmp(date, '20190528')
            fcChn  = [1:7,9:14,16];
            pulChn = [1:7,10:12,14:16] + 16;
            ppcChn = [1:7,12:13,15:16] + 32;
            vcChn  = [1:7,15:16] + 48;         
        elseif strcmp(date, '20190529')
            fcChn  = [1:7,9,11:16];
            pulChn = [1:7,10:12,14:16] + 16;
            ppcChn = [1:13,15] + 32;
            vcChn  = [1:4,6:7,15:16] + 48; 
        elseif strcmp(date, '20190530')
            fcChn  = [1:7,9:16];
            pulChn = [1:7,10:12,14:16] + 16;
            ppcChn = [1:13] + 32;
            vcChn  = [1:7,15:16] + 48;
        elseif strcmp(date, '20190531')
            fcChn  = [1:7,9:11,13:16];
            pulChn = [1:4,6,10:11,14:16] + 16;
            ppcChn = [1:13] + 32;
            vcChn  = [1:7,15:16] + 48;
        elseif strcmp(date, '20190603')
            fcChn  = [1:5,7,9:16];
            pulChn = [1:6,10:16] + 16;
            ppcChn = [1:5,7:9,11:13,15:16] + 32;
            vcChn  = [1:4,6:8,15:16] + 48; 
        elseif strcmp(date, '20190604')
            fcChn  = [1:7,9:10,12:16];
            pulChn = [1:7,10:16] + 16;
            ppcChn = [1:9,11:13,15] + 32;
            vcChn  = [1:4,6:7,15:16] + 48;             
        elseif strcmp(date, '20190607')
            fcChn  = [1:7,9:16];
            pulChn = [1:4,6:8,10:12,14:16] + 16;
            ppcChn = [1:9,11:13,15:16] + 32;
            vcChn  = [1:4,6:7,15:16] + 48;  
        elseif strcmp(date, '20190610')
            fcChn  = [1:7,9:16];
            pulChn = [1:7,10:11,14:16] + 16;
            ppcChn = [1:9,11:13,15:16] + 32;
            vcChn  = [1:10,15:16] + 48;              
            
            
        end
        validChn = {fcChn; pulChn; ppcChn; vcChn};
%%
    case '0172'
        fcChn  = 1:16;
        pulChn = [1:9,11:12,16] + 16; %excluding (17,21,22)+
        ppcChn = [1:16] + 32;
        vcChn  = [1:16] + 48;
        if  strcmp(date, '20181218') || strcmp(date, '20181219') || strcmp(date, '20181220') || strcmp(date, '20181223') || strcmp(date, '20181224')%
            pulChn = [1:9,11:13,16] + 16;
        elseif strcmp(date, '20181225') %
            fcChn  = [1:6,9:16]; %7,8 channels become noisy
            ppcChn = [1:2,4:7,9:16] + 32;
        elseif strcmp(date, '20181226') || strcmp(date, '20181228') %
            % same as default            
        elseif strcmp(date, '20181231') %
            fcChn  = [1:15];
            pulChn = [1:9,11:12] + 16;
            ppcChn = [1:7,9:16] + 32;
        elseif strcmp(date, '20190101') %
            fcChn  = [1:14];
            pulChn = [1:9,11:12] + 16;
            ppcChn = [1,3:7,10:16] + 32;
        elseif strcmp(date, '20190102') %
            fcChn  = [1:6,8:11,13:14];
            pulChn = [1:9] + 16;
            ppcChn = [3:5,9:14] + 32;
        elseif strcmp(date, '20190103') %
            fcChn  = [1:14];
            pulChn = [1:12,15:16] + 16;
            ppcChn = [1:5,8:16] + 32;    
        elseif strcmp(date, '20190104') %
            fcChn  = [1:11,13:14];
            pulChn = [1:12,15:16] + 16;
            ppcChn = [1:6,8,10:15] + 32;
        elseif strcmp(date, '20190105') %  
            pulChn = [1:12,15:16] + 16;
            ppcChn = [1,3:8,10:16] + 32;
        elseif strcmp(date, '20190107') %
            fcChn  = [1:11,13:16];
            pulChn = [1:12,15:16] + 16;
            ppcChn = [2:16] + 32; 
        elseif strcmp(date, '20190110') %
            pulChn = [1:12,15:16] + 16;
            ppcChn = [1:7,9:16] + 32;            
        elseif strcmp(date, '20190111') %
            pulChn = [1:12,15:16] + 16;
            ppcChn = [1:8,10:16] + 32; 
            vcChn  = [1:5,7,9,11:12,14,16] + 48;
        elseif strcmp(date, '20190112') %
            pulChn = [1:12,15:16] + 16;
            ppcChn = [3:8,10:16] + 32; 
            vcChn  = [1:5,7,9,11:12,14,16] + 48;
        elseif strcmp(date, '20190113') %
            pulChn = [1:12,15:16] + 16;
            ppcChn = [2:16] + 32;
            vcChn  = [1:5,8,11:12,14,16] + 48;
        elseif strcmp(date, '20190114') || strcmp(date, '20190115') || strcmp(date, '20190117')
            vcChn  = [1:5,11:12,14,16] + 48;
        elseif strcmp(date, '20190118') || strcmp(date, '20190121') 
            %same as default
        elseif strcmp(date, '20190122') 
            fcChn  = [1,2:16];
            pulChn = [1:9,11:12,15:16] + 16;
        elseif strcmp(date, '20190125') % last session of LateralVideo
            pulChn = [1:9,11:12,15:16] + 16;            
        elseif strcmp(date, '20190204') % start MozartAudio
            %same as default
        elseif strcmp(date, '20190206')
            pulChn = [1:12] + 16;
        elseif strcmp(date, '20190208')
            fcChn  = [1:10,11:16];
            pulChn = [1:12] + 16;
        elseif strcmp(date, '20190212')
            fcChn  = [1:12,14:16];
            pulChn = [1:9,11:12] + 16;
            ppcChn = [1,14:16] + 32;
            vcChn  = [4,11:16] + 48;
        elseif strcmp(date, '20190220')
            fcChn  = [1:12,14:16];
            pulChn = [1:9,11:12] + 16;
            ppcChn = [1:13] + 32;
        end
        validChn = {fcChn; pulChn; ppcChn; vcChn};
    case '0168'
        if strcmp(date, '20180629') % recording 003 on 6/29
            vcChn  = 1:16;
            ppcChn = sort(chanOrder32(1,1:15)) + 16;
            pulChn = [1:16] + 16 + 32;
        elseif datetime(date, 'InputFormat', 'yyyyMMdd') <= datetime('20180712','InputFormat', 'yyyyMMdd') % half PPC channels are faulty
            pulChn = 1:16;
            ppcChn = sort(chanOrder32(1,1:15)) + 16;
            vcChn  = [1:16] + 16 + 32;
        elseif datetime(date, 'InputFormat', 'yyyyMMdd') <= datetime('20180807','InputFormat', 'yyyyMMdd') % half PPC channels are faulty
            pulChn = 1:16;
            ppcChn = [1:8, 11:22, 24:32] + 16;% 7/13 switch to a better PPC headstage
            vcChn  = [1:16] + 16 + 32;
        elseif datetime(date, 'InputFormat', 'yyyyMMdd') >= datetime('20180808','InputFormat', 'yyyyMMdd') % half PPC channels are faulty
            pulChn = 1:16; %according to LFP
            ppcChn = [1:7,25:32] + 16;
            vcChn  = [1:16] + 16 + 32;            
        else
            pulChn = [1, 3:16];
            ppcChn = [1:8, 12:22, 24:32] + 16;
            vcChn  = [1:16] + 16 + 32;

        end
        validChn = {pulChn; ppcChn; vcChn}; %original channel order from 1 to 64 but only valid ones

    case '0169'
        lpulChn   = [1:2,6:7,9:12,14:16];
        lppcChn   = [1:16] + 16 ;
        rppcChn   = [1:11,14:16] + 16 + 16;
        rpulChn   = [3:6,8:9,11:16] + 16 + 16 + 16;         
        lvcEEGChn = [1]  + 16 + 16 + 16 + 16;
        rvcEEGChn = [1]  + 16 + 16 + 16 + 16 + 1;
        switch date
            case '20180713'
                lpulChn = [1:3,10:12,14:16]; % exclude alpha channels
                lppcChn = [1:4,6:7,9:16] + 16;
                rpulChn = [2:16] + 16 + 16 + 16; 
            case '20180716'
                lpulChn = [1:3,10:12,14:16];
                lppcChn = [1:5,8:16] + 16;
                rppcChn = [1:10,12:16] + 16 + 16;
                rpulChn = [1:3,6:12,14:15] + 16 + 16 + 16;
            case '20180717'
                lpulChn = [1:3,7:12,14:16];
                lppcChn = [2:6,8:16] + 16;
                rppcChn = [2:6,8:11,14:16] + 16 + 16;
                rpulChn = [1:3,6:12,14:15] + 16 + 16 + 16; 
            case '20180718'
                lpulChn = [1:2,7:16];
                lppcChn = [1:2,4:6,8:16] + 16;
                rppcChn = [2:3,5:9,12:15] + 16 + 16;
                rpulChn = [1:2,5:12,14:15] + 16 + 16 + 16; 
            case '20180719'
                lpulChn = [1:2,7:16];
                lppcChn = [2:6,8:10,12:16] + 16;
                rppcChn = [3:8,12:15] + 16 + 16;
                rpulChn = [1:2,5:12,14] + 16 + 16 + 16; 
            case '20180720'
                lpulChn = [1:2,7:16];
                lppcChn = [1:2,4:10,12:16] + 16;
                rppcChn = [1:7,13:15] + 16 + 16;
                rpulChn = [1:3,5:11,15] + 16 + 16 + 16;
            case '20180723'    
                lpulChn = [1:2,7:16];
                lppcChn = [1:5, 8:13,15:16] + 16;
                rppcChn = [2:7, 13:14] + 16 + 16;
                rpulChn = [1:2,6:12,14] + 16 + 16 + 16;
            case '20180724'
                lpulChn = [1:2,6:16];
                lppcChn = [1,4:6, 8:16] + 16;
                rppcChn = [1:14] + 16 + 16;
                rpulChn = [1:2,6:12,14] + 16 + 16 + 16;
            case '20180725' %023 is bilateral alpha session
                lpulChn = [1:2,6:16];
                lppcChn = [4,7:16] + 16;%[9,11:16] + 16; %PPC like spectrogram
                rppcChn = [1:9,16] + 16 + 16;
                rpulChn = [1:2,6:12,14] + 16 + 16 + 16;
            case '20180726' %full video look good
                lpulChn = [1:2,6:7,9:14,16];
                lppcChn = [1:6, 8:16] + 16;
                rppcChn = [1:10, 14:16] + 16 + 16;
                rpulChn = [1:4,9:14] + 16 + 16 + 16;
            case '20180727'
                lpulChn = [1:2,6:7,9:14,16];
                lppcChn = [1:6, 8:16] + 16;
                rppcChn = [1:10, 14:16] + 16 + 16;
                rpulChn = [1:5, 13:16] + 16 + 16 + 16;
            case '20180807'
                lpulChn = [1:2,6:7,9:14,16];
                lppcChn = [1:6, 8:16] + 16;
                rppcChn = [1:10, 14:16] + 16 + 16;
                rpulChn = [1:5, 13:16] + 16 + 16 + 16; 
            case '20180808'
                lpulChn   = [1, 6:14, 16];
                lppcChn   = [2:6,8:16] + 16 ;
                rpulChn   = [5:6,10:11,13] + 16 + 16 + 16;
            case '20180809'
                lpulChn   = [1:2,6:14,16];
                rppcChn   = [1:11,14:16] + 16 + 16;
                rpulChn   = [2:4,13,16] + 16 + 16 + 16;         
            case '20180814'
                lpulChn   = [1:7,9:11,15:16];
                rppcChn   = [1:11,14:16] + 16 + 16;
                rpulChn   = [1:10,14:15] + 16 + 16 + 16;        
            case '20180815' 
                lpulChn   = [1:14,16];
                rppcChn   = [1:4,6:16] + 16 + 16;
                rpulChn   = [1:9,11:16] + 16 + 16 + 16;            
            case '20180816' %full view session
                lpulChn   = [1:2,6:7,9:12,14:16]; %show LPl-like theta
                rppcChn   = [1:11,14:16] + 16 + 16;
                rpulChn   = [3:6,8:9,11:16] + 16 + 16 + 16;
            case '20180820' %full view session
                lpulChn   = [1,6:7,9:12,14:16];
                rppcChn   = [1:11,14:16] + 16 + 16;
                rpulChn   = [2:7,13:14] + 16 + 16 + 16;
            case '20180822' %full view session checked lfp PPC effect very good
                lpulChn   = [1,6:7,9:12,14:16];
                rppcChn   = [1:2,4:11,14:16] + 16 + 16;
                rpulChn   = [1:4,6:8,10:14,16] + 16 + 16 + 16;
            case '20180831'
                lpulChn   = [1,6:7,9:12,14:16];
                rppcChn   = [1:11,14:16] + 16 + 16;
                rpulChn   = [3:5,8:9,11:16] + 16 + 16 + 16;
            case '20180904'
                lpulChn   = [1:2,6:7,9:12,15:16];
                rppcChn   = [1:11,14:16] + 16 + 16;
                rpulChn   = [3:6,10:16] + 16 + 16 + 16; 
            case '20180905'
                lpulChn   = [7,9:12,15:16]; %090 session
                lppcChn   = [9:16] + 16;
                rppcChn   = [8:11,14:16] + 16 + 16;
                rpulChn   = [1:6,14:16] + 16 + 16 + 16; 
            case '20180906'
                lpulChn   = [1,6:7,9:12,16];
                rppcChn   = [1:11,14:15] + 16 + 16;
                rpulChn   = [3:6,9:11,14:16] + 16 + 16 + 16; 
            case '20180910'
                lpulChn   = [1:2,6:7,9:11,14:16];
                rppcChn   = [1:11,14:16] + 16 + 16;
                rpulChn   = [3:6,8:9,11:16] + 16 + 16 + 16;
        end
        
        validChn  = {lpulChn; lppcChn; rppcChn; rpulChn; lvcEEGChn; rvcEEGChn}; %original channel order from 1 to 64 but only valid ones

end

reorderedChn = {};
numPrevChn   = 0;
for iRegion = 1:length(validChn)            
    reorderedChn{iRegion} = [1:length(validChn{iRegion})]+numPrevChn;
    numPrevChn = numPrevChn + length(validChn{iRegion});
end

end


