function scaleFigure(hFig, yscale, xscale)


pos = get(hFig, 'Position'); % [left bottom width height]

set(hFig, 'Position', pos.*[1, 1, xscale, yscale]);

end

