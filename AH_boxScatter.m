function h = AH_boxScatter(data,groupID,xLabel,setOrder,alpha,plotLine)
% This funtion takes in arrays of data and groupID, then plot boxplot and
% scatter plot on top
% AH 2/13/2020

if nargin <= 3
    setOrder = 'sorted';
    alpha = 0.4;
    plotLine = 0;
elseif nargin <= 4
    alpha = 0.4;
    plotLine = 0;
end

% Needs to run scatter -> line -> box, if run box first, it will generate 5 diff regression line fit
[C, ~, ic]= unique(groupID,setOrder); % sorted order
xratio = max(ic)/max(groupID); % ratio of new x over old x, plot will appear 0-1 but internally box plot see it as 1-5
scatter(ic,data,'filled','MarkerFaceAlpha',alpha','jitter','on','jitterAmount',0.15);
%scatter(groupID*5,data,'filled','MarkerFaceAlpha',alpha','jitter','on','jitterAmount',0.03);

%% Draw best fitted line (linear least squares fit) -- OLD, use bend-corr line below
% Method 1 (if using ic,data): 
%h.line = lsline;
% Add this if using groupID, data):
%margin = 0.1; % add left/right margin for boxplot
%xlim([min(groupID)-margin,max(groupID)+margin]);

hold on;

% Pearson
[h.rPearson h.pPearson] = corrcoef(groupID,data,'Rows','pairwise'); 
h.rPtext = sprintf('%0.3f', h.rPearson(1,2)); h.pPtext = sprintf('%0.3f', h.pPearson(1,2));
% Spearman
[h.rSpearman,h.pSpearman] = corr(groupID,data,'Rows','pairwise','Type','Spearman'); 
h.rStext = sprintf('%0.3f',h.rSpearman); h.pStext = sprintf('%0.3f', h.pSpearman);
% Bend-correlation
[h.rBend,t,h.pBend,~,~,~,~,line] = bendcorr(groupID,data,0); % 0 to suppress plotting
h.rBtext = sprintf('%0.3f',h.rBend); h.pBtext = sprintf('%0.3f', h.pBend);

if plotLine
    % Plot bend-correlation fitted line
    h.Bline = refline(line.slope/xratio,line.intercept); 
    h.BlineMeta.slope = line.slope;
    h.BlineMeta.intercept = line.intercept;
    set(h.Bline,'Color','r','LineWidth',2);

    % Plot CI and fill
    y1 = refline(line.CIslope(1)/xratio,line.CIintercept(1)); set(y1,'Color','r');
    y2 = refline(line.CIslope(2)/xratio,line.CIintercept(2)); set(y2,'Color','r');
    y1 = get(y1); y2 = get(y2);
    xpoints=[[y1.XData(1):y1.XData(2)],[y2.XData(2):-1:y2.XData(1)]];
    step1 = y1.YData(2)-y1.YData(1); step1 = step1 / (y1.XData(2)-y1.XData(1));
    step2 = y2.YData(2)-y2.YData(1); step2 = step2 / (y2.XData(2)-y2.XData(1));
    filled=[[y1.YData(1):step1:y1.YData(2)],[y2.YData(2):-step2:y2.YData(1)]];
    hold on; fillhandle=fill(xpoints,filled,[1 0 0]);
    set(fillhandle,'EdgeColor',[1 0 0],'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
    h.BlineMeta.CIslope = line.CIslope;
    h.BlineMeta.CIintercept = line.CIintercept;
end

h.titleText = {['Pearson r=' h.rPtext ' p=' h.pPtext];...
         ['Spearman r=' h.rStext ' p=' h.pStext];...
         ['Bend-corr r=' h.rBtext ' p=' h.pBtext]}; 
     
% Plot box the last
h.box = boxplot(data,groupID,'Labels',xLabel,'Notch','on','medianstyle','line'); 
% To hide outliers
%set(h.box(7,:),'Visible','off');

end