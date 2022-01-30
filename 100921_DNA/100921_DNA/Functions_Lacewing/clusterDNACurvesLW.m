function [coord_cell,coord_vect,K] = clusterDNACurvesLW(in,pixel_on,K)
%findRegionLW performs a DBSCAN clustering on the ISFET array to extract
% sensing regions.
% Inputs :
% - in: the input vector for clustering.
% - pixel_on : vector for active pixels used to reshape in into the array.

% First convert array into row format
data_3 = array_to_rowLW(in,pixel_on);
rng(1);

% Determine number of clusters
sse = zeros(1,K);
for k = 1:K
    [~,~,dist] = kmeans(data_3,k);
    sse(k) = sum(dist);
end
sse = sse/sse(1);
x = (1:K)/K;
% figclf(10); plot(x,sse); hold all;
pol = polyfit(x([1,K]),sse([1,K]),1);
y_pol = polyval(pol,x); 
%     plot(x,y_pol); axis equal;
dist_line = zeros(K,1);
for k = 1:K
    pol_perp = [-1/pol(1) sse(k)+1/pol(1)*x(k)];
%         y_pol_perp = polyval(pol_perp,x);
%         plot(x,y_pol_perp);
    % Point on line
    pt1_x = (pol_perp(2)-pol(2))/(pol(1)-pol_perp(1));
    pt1_y = polyval(pol,pt1_x);
    dist_line(k) = norm([pt1_x,pt1_y] - [x(k),sse(k)]);
end

[~,K] = max(dist_line);
coord_vect = kmeans(data_3,K);
for i = 1:K
    coord_cell{i} = find(coord_vect==i);
end