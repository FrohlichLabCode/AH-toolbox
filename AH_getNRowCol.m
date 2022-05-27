function [nRow,nCol] = AH_getNRowCol(fig)
% This function is to get number of rows and cols from a regular grid
% figure with handle fig
N=numel(fig.Children);
counter = 0;
for n = 2:2:N % skip colorbar
    counter = counter+1;
    pos1(counter) = fig.Children(n).Position(1);
    pos2(counter) = fig.Children(n).Position(2);
end
nCol = numel(unique(pos1));
nRow = numel(unique(pos2));
end