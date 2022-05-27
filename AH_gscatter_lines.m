function h = AH_gscatter_lines(x,y,thisLabel,labelTypeIDs,color1)
% This script will plot scatter plot, colored by group
% And calulate correlation coeffient for entire data and each group
% seperately, and then comparing 2 groups.

gscatter(x,y, thisLabel,color1,'os',4); % x,y, size, color % can't add transparency to gscatter https://www.mathworks.com/matlabcentral/answers/574957-how-can-i-make-the-plot-transparent-in-a-gscattter?s_tid=srchtitle
% add jitter -- doesn't work
%gscatter(x,y, thisLabel,color1,'os',4, 'jitter','on', 'jitterAmount',0.5); % x,y, size, color % can't add transparency to gscatter https://www.mathworks.com/matlabcentral/answers/574957-how-can-i-make-the-plot-transparent-in-a-gscattter?s_tid=srchtitle
h.hl = lsline;
% Compare correlations
h.P=[];
h.S=[];
h.B=[];
h.rPtestAll=[];
h.pPtestAll=[];
h.rStestAll=[];
h.pStestAll=[];
h.rBtestAll=[];
h.pBtestAll=[];

for i = 1:numel(labelTypeIDs) % for each line, get slope (correlation)
    if strcmp(labelTypeIDs(i), 0)
        continue
    end
    % Get mask for each label (eg. opto condition)
    mask(:,i) = thisLabel == labelTypeIDs(i);
    [h.P(i,1) h.P(i,2)] = corr(x(mask(:,i)),y(mask(:,i)),'Rows','pairwise','Type','Pearson');
    [h.S(i,1) h.S(i,2)] = corr(x(mask(:,i)),y(mask(:,i)),'Rows','pairwise','Type','Spearman');
    [h.B(i,1),~,h.B(i,2),~,~,~,~,h.Bline{i}] = bendcorr(x(mask(:,i)),y(mask(:,i)),0);
    h.rPtest{i} = sprintf('%0.2f', h.P(i,1)); h.pPtest{i} = sprintf('%0.2f', h.P(i,2));
    h.rPtestAll = [h.rPtestAll ' ' h.rPtest{i}]; h.pPtestAll = [h.pPtestAll ' ' h.pPtest{i}];    
    h.rStest{i} = sprintf('%0.2f', h.S(i,1)); h.pStest{i} = sprintf('%0.2f', h.S(i,2));
    h.rStestAll = [h.rStestAll ' ' h.rStest{i}]; h.pStestAll = [h.pStestAll ' ' h.pStest{i}];
    h.rBtest{i} = sprintf('%0.2f', h.B(i,1)); h.pBtest{i} = sprintf('%0.2f', h.B(i,2));
    h.rBtestAll = [h.rBtestAll ' ' h.rBtest{i}]; h.pBtestAll = [h.pBtestAll ' ' h.pBtest{i}];
    
end

[h.p12, z, za, zb] = corr_rtest(h.B(1,1), h.B(2,1), sum(mask(:,1)), sum(mask(:,2))); h.p12text = sprintf('%0.2f', h.p12(2)); % p(2) is pvalue of two-tailed test

% Pearson
[h.rPearson h.pPearson] = corrcoef(x,y,'Rows','pairwise'); 
h.rPtext = sprintf('%0.2f', h.rPearson(1,2)); h.pPtext = sprintf('%0.2f', h.pPearson(1,2));
% Spearman
[h.rSpearman,h.pSpearman] = corr(x,y,'Rows','pairwise','Type','Spearman'); 
h.rStext = sprintf('%0.2f',h.rSpearman); h.pStext = sprintf('%0.2f', h.pSpearman);
% Bend-correlation (don't need CI)
[h.rBend,~,h.pBend] = bendcorr(x,y,0); % 0 to suppress plotting
h.rBtext = sprintf('%0.2f',h.rBend); h.pBtext = sprintf('%0.2f', h.pBend);

h.titleText = {['All rP=' h.rPtext ' pP=' h.pPtext ', rS=' h.rStext ' pS=' h.pStext ', rB=' h.rBtext ' pB=' h.pBtext];...
    ['rP=' h.rPtestAll ' pP=' h.pPtestAll ' p12=' h.p12text];... 
    ['rS=' h.rStestAll ' pS=' h.pStestAll ' p12=' h.p12text];... 
    ['rB=' h.rBtestAll ' pB=' h.pBtestAll ' p12=' h.p12text]}; 
            