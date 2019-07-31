function fig = AH_figure(numRows, numCols, name)
%fig = AH_figure(numRows, numCols, name); %numRows, numCols, name

%screensize = get( groot, 'Screensize' );
if ~exist('name','var'); name = 'figure'; end
fig = figure('name',name,'Position',[10 50 320*numRows 270*numCols]);%x,y,width,height
end
