function running_intersect = mintersect(varargin)
% Multiple set intersection.

    if isempty(varargin)
        error('Not enough inputs specified.')
    end

    running_intersect = varargin{1};
    for i = 2:nargin
        running_intersect = intersect(varargin{i}, running_intersect);
    end

end