function heatmap_cluster(data, rownames, colnames, datarange, cmap)
	% Create heatmap based on matrix stored in data
	% Input arguments: data - data matrix
	% rownames - names of rows
	% colnames - names of columns
	% datarange - range of values that will be associated with a color
	% cmap - colormap used to create the plot
    data(isnan(data)) = -1461;
    D = pdist(data);
    tree = linkage(D, 'average');
    leafOrder_row = optimalleaforder(tree, D);
    D = pdist(data');
    tree = linkage(D, 'average');
    leafOrder_column = optimalleaforder(tree, D);
	data(data == -1461) = NaN;
	
	if nargin < 5
		cmap = colormap;
	end

    if length(rownames) >= 1
        heatmap(data(leafOrder_row, leafOrder_column), colnames(leafOrder_column), ...
            rownames(leafOrder_row), [], 'GridLine', ':', 'ShowAllTicks', true, 'TickAngle', 45, ...
			'MinColorValue', datarange(1), 'MaxColorValue', datarange(2), 'NaNColor', [1 1 1],...
			'Colormap', cmap);
    else
        heatmap(data(leafOrder_row, leafOrder_column), [], [], [], 'GridLine', ':', ...
            'ShowAllTicks', true, 'TickAngle', 45, 'MinColorValue', datarange(1), 'MaxColorValue', ...
            datarange(2), 'NaNColor', [1 1 1],'Colormap', cmap);
    end

end
