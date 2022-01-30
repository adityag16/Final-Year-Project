function [coord_cluster_final,res_cluster_final,tm_cluster_final, ...
    density_cluster_final, res_cluster, tm_cluster, density_cluster] ...
    = classifyClusterLW(coord_cluster,time,pixel_on,res_exp_log,tm,x_stable)
%classifyClusterLW extracts relevant metrics from a set of clusters on the
% TTN array and classifies according to 'active' and 'barrier' clusters

% Constants
[ROWS,COLS] = getConstantsLW;

% Initialisation of metrics
C = length(coord_cluster);
n_cluster = cellfun(@length,coord_cluster)';
res_cluster = zeros(C,1); % Residual from exp approximation
tm_cluster = zeros(C,1); % Average peak derivative time
density_cluster = zeros(C,1); % Cluster density within a window
res_cluster_final = zeros(2,1); % Residual from exp approximation
tm_cluster_final = zeros(2,1); % Average peak derivative time
density_cluster_final = zeros(2,1); % Cluster density within a window

% Metrics for each cluster
for c = 1:C
    res_cluster(c) = mean(res_exp_log(pixel_on(coord_cluster{c})));
    tm_cluster(c) = mean(tm(pixel_on(coord_cluster{c})))/60;
end
% Sub-window for density calculation
array_wdo_x = (15:COLS)';
array_wdo_y = (15:50)';
Wx = length(array_wdo_x);
Wy = length(array_wdo_y);
array_wdo = [repelem(array_wdo_x,Wy,1),repmat(array_wdo_y,Wx,1)];
array_wdo_v = a2v(array_wdo,ROWS,COLS);
array_wdo_on = intersect(array_wdo_v,pixel_on);
Won = length(array_wdo_on);
% Density of each cluster
for c = 1:C
    idx = intersect(array_wdo_on,pixel_on(coord_cluster{c}));
    density_cluster(c) = length(idx)/Won;
end
    
% Select clusters based on all metrics
idx_res = find(res_cluster>-1000);
idx_density = find(density_cluster>.2);
idx_tm = find(tm_cluster>time(x_stable)/60+1);
% Criteria for cluster classification
idx_cluster = intersect(idx_res,intersect(idx_tm,idx_density));
idx_cluster_n = setdiff((1:C)',idx_cluster);
coord_cluster_final{1} = cat(1,coord_cluster{idx_cluster}); % On cluster
coord_cluster_final{2} = cat(1,coord_cluster{idx_cluster_n}); % Barrier cluster

% Final metric vectors
res_cluster_final(1) = sum(n_cluster(idx_cluster).*res_cluster(idx_cluster))./sum(n_cluster(idx_cluster));
res_cluster_final(2) = sum(n_cluster(idx_cluster_n).*res_cluster(idx_cluster_n))/sum(n_cluster(idx_cluster_n));
tm_cluster_final(1) = sum(n_cluster(idx_cluster).*tm_cluster(idx_cluster))./sum(n_cluster(idx_cluster));
tm_cluster_final(2) = sum(n_cluster(idx_cluster_n).*tm_cluster(idx_cluster_n))/sum(n_cluster(idx_cluster_n));
density_cluster_final(1) = sum(density_cluster(idx_cluster));
density_cluster_final(2) = sum(density_cluster(idx_cluster_n));
end

