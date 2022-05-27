function hBar = AH_plotTableAsGroupedBar(barTable, colLabel, displayDigit, errTable, data, masks, rowLabel) % Table, xTickLabel, displayDigit
% Note that the order of the rows in table should match the xTickLabel
% From errTable onwards are optional variables.
% AH 2/7/2020
% updated version

if ~exist('errTable','var')
    errTable = [];
end
% data is optional, if want to add each data point
if ~exist('data','var') 
    data = [];
end
if ~exist('masks','var') 
    masks = [];
end
if ~exist('rowLabel','var')
    rowLabel=[];    
end
    [nRow, nCol] = size(barTable); % nCol=number of big group, nRow=bars within a group
    % Plot Bar graph
    if size(barTable,1) == 1 % only 1 row, not grouped, then plot each bar at a time to get diff color
        for iCol = 1:nCol
            hBar(1,iCol) = bar(iCol, table2array(barTable(:,iCol)));
            %set(gca,'XTickLabel',xTickLabel{iCol}); % doesn't work, need to put label for each column
            hold on
        end
    elseif size(barTable,2) == 1 % only 1 row, not grouped, then plot each bar at a time to get diff color
        for iRow = 1:nRow
            hBar(iRow,1) = bar(iRow, table2array(barTable(iRow,:)));
            %set(gca,'XTickLabel',xTickLabel{iCol}); % doesn't work, need to put label for each column
            hold on
        end
        
    else % Grouped by column, each column is one color automatically
        hBar = bar(table2array(barTable));
    end
    % don't use hBar = bar(categorical(acc{1:4,2}), table2array(acc(1:4,trialTypes))); it will rearrange the group order which doesn't match the name

    set(gca,'XTickLabel',colLabel);
    for iCol = 1:nCol % for each column of graph, add text of value to it
        if ~isempty(displayDigit) % if empty, don't display digit
            if displayDigit == 0
                text([1:1:nRow]+0.23*(iCol-2),zeros(1,nRow),string(round(barTable{:,iCol})'),'horizontalalignment','center','verticalalignment','bottom','FontSize', 8)
            elseif displayDigit > 0
                if nCol == 1
                    text([1:1:nRow],zeros(1,nRow),string(round(barTable{:,iCol},displayDigit)'),'horizontalalignment','center','verticalalignment','bottom','FontSize', 8) 
                elseif nCol == 4
                    text([1:1:nRow]+0.17*(iCol-2.5),zeros(1,nRow),string(round(barTable{:,iCol},displayDigit)'),'horizontalalignment','center','verticalalignment','bottom','FontSize', 8)                    
                elseif mod(nCol,2) == 1 % odd number of rows, center at 0
                    text([1:1:nRow]+0.23*(iCol-2),[-12].*ones(1,nRow),string(round(barTable{:,iCol},displayDigit)'),'horizontalalignment','center','verticalalignment','top','FontSize', 8)
                else % even number
                    if nRow == 1
                        text(iCol,0,string(round(barTable{:,iCol},displayDigit)'),'horizontalalignment','center','verticalalignment','top','FontSize', 8)
                    else
                        text([1:1:nRow]+0.3*(iCol-1.5),[-12].*ones(1,nRow),string(round(barTable{:,iCol},displayDigit)'),'horizontalalignment','center','verticalalignment','top','FontSize', 8)
                    end
                end
            end
        end
    end
    
    
    hold on
    if ~isempty(errTable) % add error bar
        for iCol = 1:nCol 
            neg = errTable{:,iCol}.*(barTable{:,iCol}<0); % plot bottom half if mean bar is negative
            pos = errTable{:,iCol}.*(barTable{:,iCol}>=0); 
            if nCol == 1
                errorbar([1:1:nRow],barTable{:,iCol},neg,pos,'k.','CapSize',round(24/nCol))
            elseif nCol == 4 % first number is the spread of sub errorbars, 2nd number is spread of main groups
                errorbar([1:1:nRow]+0.18*(iCol-2.5),barTable{:,iCol},neg,pos,'k.')
            elseif mod(nCol,2) == 1
                errorbar([1:1:nRow]+0.22*(iCol-2),barTable{:,iCol},neg,pos,'k.')
            else
                if nRow == 1
                    errorbar(iCol,barTable{:,iCol},neg,pos,'k.','CapSize',round(48/nCol)) % 22
%                 elseif nCol == 1
%                     errorbar([1:1:nRow]+0.29*(iCol-1.5),barTable{:,iCol},neg,pos,'k.','CapSize',18) % 22
                else
                    errorbar([1:1:nRow]+0.29*(iCol-1.5),barTable{:,iCol},neg,pos,'k.')
                end
            end
        end
    end
    
    if ~isempty(data) % add data points
        if nRow == 1
            nColor = max(nCol,size(hBar,2));  % was size(hBar,2)
        else
            nColor = nRow;
        end
        for iHM = 1:nRow
            if istable(data)==1
                mask = masks(:,iHM);
                nSess = sum(mask);
            else % don't need mask
                nSess = size(data,1)/nCol; 
                %nSess = size(data,1); % for SUPSTH
            end
            if nRow == 1
                x = reshape(repmat([1:nColor],[nSess,1]),1,[]);
            else
                if nCol == 3
                    x = reshape(repmat([iHM-0.22,iHM,iHM+0.22],[nSess,1]),1,[]);
                elseif nCol == 2
                    x = reshape(repmat([iHM-0.15,iHM+0.15],[nSess,1]),1,[]);
                elseif nCol == 1                
                    x = reshape(repmat([iHM],[nSess,1]),1,[]);
                elseif nCol == 4 % off set from the center of each subgroup center
                    x = reshape(repmat([iHM-0.28,iHM-0.12,iHM+0.12,iHM+0.28],[nSess,1]),1,[]);
                end
            end
            if istable(data)==1 % check if data is a table
                y = reshape(data{mask,rowLabel},1,[]);
            else % data is one long vector, use mask to select which ones to include
                y = data(:,iHM); % note data has to be stacked (11112222) to match x
            end
            sz = 15;
            shade = 1/2; % make points darker, smaller=darker
            %shade = 1; %original color, blend in with bar
            c = []; % form nx3 rgb color vector
            if nCol == 1
                c = repmat(hBar(iHM).FaceColor*shade,[nSess,1]);
            else
                for iColor = 1:nCol % or nColor?
                    c = [c; repmat(hBar(iColor).FaceColor*shade,[nSess,1])];
                end
            end
            %group = reshape(repmat([1,2,3],[nSess,1]),1,[]); % gscatter doesn't have jitter
            if nRow == 1 || nCol == 1
                jitterAmount = 0.2;
            else
                jitterAmount = 0.05;
            end
            scatter(x,y,sz,c,'filled','MarkerFaceAlpha',0.3','jitter','on','jitterAmount',jitterAmount);
        end 
    end
end