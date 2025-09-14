function set_background(x,y,colors,thresholds)
idx = ones(1,numel(x));
thresholds = [-Inf thresholds Inf];
for ii = 2:numel(thresholds)
    idx(y >= thresholds(ii-1) & y < thresholds(ii)) = ii-1;
end
yl = get(gca(),'YLim').';
x = (x(1:end-1)+x(2:end))/2;
x = [2*x(1)-x(2) x 2*x(end)-x(end-1)];
p = patch( ...
    'XData',x((2:end)+[0;0;-1;-1;0]), ...
    'YData',repmat(yl([1 2 2 1 1]),1,numel(x)-1), ...
    'FaceColor','flat', ...
    'FaceVertexCData',colors(idx,:), ...
    'EdgeColor','none');
ch = get(gca(),'Children');
ch(ch == p) = [];
set(gca(),'Children',[ch(:); p],'Layer','top');
end