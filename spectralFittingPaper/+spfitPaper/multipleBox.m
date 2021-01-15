function ax = multipleBox(X,grouping,varargin)
% multipleBox  Produce multiple box plots stacked on top of each other.
%   multipleBox(X,grouping,varargin) produces multiple box plots, looping
%   over the third dimension. The first two are passed to boxplot as
%   normal.
%
%   See also boxplot.
rows = size(X,3);

BigAx = newplot(gca);
pos = get(BigAx,'Position');
height = pos(4)/rows;
space = 0.02;
ax = gobjects(rows,1);

for jj = rows:-1:1
    ax(jj) = subplot(rows,1,jj);
    boxplot(squeeze(X(:,:,jj)),grouping,varargin{:})
    if jj ~= rows
        set(gca, 'XTickLabel', []);
    end
    ax(jj).Position = [pos(1) pos(2)+(rows-jj)*height ...
            pos(3) height*(1-space)];
end

end