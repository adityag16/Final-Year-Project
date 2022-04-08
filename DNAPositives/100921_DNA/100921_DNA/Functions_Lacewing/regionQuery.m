%------------------------------------------------------------------------
% Find indices of all points within distance E of point with index n
% This function could make use of a precomputed distance table to avoid
% repeated distance calculations, however this would require N^2 storage.
% Not a big problem either way if the number of points being clustered is
% small.   For large datasets this function will need to be optimised.

% Arguments:
%              P - the dim x Npts array of data points
%              n - Index of point of interest
%              E - Distance threshold

function neighbours = regionQuery(P, n, E)
    [dim, Npts] = size(P);
    neighbours = [];
    
    for i = 1:Npts
        if i ~= n
            % Test if distance^2 < E^2 
            v = P(:,i)-P(:,n);
            dist2 = v(1:2)'*v(1:2);
            if (dist2 < E(1)^2) && ((abs(v(3)) < abs(E(2))))
               neighbours = [neighbours i];     
            end
        end
    end
    
end % of regionQuery
