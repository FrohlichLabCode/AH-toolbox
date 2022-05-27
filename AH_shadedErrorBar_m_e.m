function l1 = AH_shadedErrorBar_m_e(x, y, dy, colorstyle)
% l1 is handle of the line
% eg. to display legend: legend([l1 l2],'theta','gamma','Location','east');
% x and y has to be column vector
% AH 2020/5/7

h = fill([x;flipud(x)],[y-dy;flipud(y+dy)],colorstyle(1),'linestyle','none');
hold on
set(h,'facealpha',.5); % make translucent
l1 = plot(x,y,colorstyle);