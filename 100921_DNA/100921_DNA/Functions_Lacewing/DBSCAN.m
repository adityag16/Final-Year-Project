% DBSCAN DBSCAN clustering algorithm
%
% Usage:  [C, ptsC, centres] = dbscan(P, E, minPts)
%
% Arguments:
%         P - dim x Npts array of points.
%         E - Distance threshold.
%    minPts - Minimum number of points required to form a cluster.
%
% Returns:
%         C - Cell array of length Nc listing indices of points associated with
%             each cluster.
%      ptsC - Array of length Npts listing the cluster number associated with
%             each point.  If a point is denoted as noise (not enough nearby
%             elements to form a cluster) its cluster number is 0.
%   centres - dim x Nc array of the average centre of each cluster.

function [C, ptsC, centres] = DBSCAN(P, E, minPts)
    
    [dim, Npts] = size(P);
    
    ptsC  = zeros(Npts,1);
    C     = {};
    Nc    = 0;               % Cluster counter.
    Pvisit = zeros(Npts,1);  % Array to keep track of points that have been visited.
    
    for n = 1:Npts
       if (~Pvisit(n))                          % If this point not visited yet
           Pvisit(n) = 1;                       % mark as visited
           neighbourPts = regionQuery(P, n, E); % and find its neighbours

           if length(neighbourPts) < minPts-1  % Not enough points to form a cluster
               ptsC(n) = 0;                    % Mark point n as noise.
           
           else                % Form a cluster...
               Nc = Nc + 1;    % Increment number of clusters and process
                               % neighbourhood.
           
               C{Nc} = [n];    % Initialise cluster Nc with point n
               ptsC(n) = Nc;   % and mark point n as being a member of cluster Nc.
               
               ind = 1;        % Initialise index into neighbourPts array.
               
               % For each point P' in neighbourPts ...
               while ind <= length(neighbourPts)
                   
                   nb = neighbourPts(ind);
                   
                   if ~Pvisit(nb)        % If this neighbour has not been visited
                       Pvisit(nb) = 1;   % mark it as visited.
                       
                       % Find the neighbours of this neighbour and if it has
                       % enough neighbours add them to the neighbourPts list
                       neighbourPtsP = regionQuery(P, nb, E);
                       if length(neighbourPtsP) >= minPts
                           neighbourPts = [neighbourPts , neighbourPtsP];
                       end
                   end            
                   
                   % If this neighbour nb not yet a member of any cluster add it
                   % to this cluster.
                   if ~ptsC(nb)  
                       C{Nc} = [C{Nc} nb];
                       ptsC(nb) = Nc;
                   end
                   
                   ind = ind + 1;  % Increment neighbour point index and process
                                   % next neighbour
               end
           end
       end
    end
    
    % Find centres of each cluster
    centres = zeros(dim,length(C));
    for n = 1:length(C)
        for k = 1:length(C{n})
            centres(:,n) = centres(:,n) + P(:,C{n}(k));
        end
        centres(:,n) = centres(:,n)/length(C{n});
    end

end % of dbscan    