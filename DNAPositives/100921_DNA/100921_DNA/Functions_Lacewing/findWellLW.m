function [coord_well] = findWellLW(clusters,idx,pixel_on)
%FINDWELLLW finds the indices of the main well from the input clusters and 
% the index of the main cluster idx.

% Constants
[ROWS,COLS,P,coord_chem,ROWS_temp,COLS_temp,npt,coord_temp,coord_peltier] = getConstantsLW;

% [row col] coordinates of the main cluster - pixels to include
coord_main = clusters{idx};

% [row col] coordinates of the other clusters - 'barrier' pixels
idx_n = setdiff(1:length(clusters),idx);
coord_barrier = vertcat(clusters{idx_n});

% Note: temp and off pixels are ignored and considered 'neutral'

% Generate array where off pixels = 1, main cluster = 1 and others = 2
array_vect=zeros(1,P);
array_vect(coord_temp) = 1;
array_vect(pixel_on(coord_main)) = 1;
array_vect(pixel_on(coord_barrier)) = 2;
array = reshape(array_vect,ROWS,COLS);

% Starting point    
start_pt = floor(median(v2a(pixel_on(coord_main),ROWS,COLS)));
while array(start_pt(1),start_pt(2))~=1
    if start_pt(1)<COLS
        start_pt(1) = start_pt(1)+1;
    else
        start_pt(2) = start_pt(2)+1;
    end
end
coord_region = findRegionLW(array,start_pt);

% Output
coord_well = intersect(coord_region,pixel_on);
end

