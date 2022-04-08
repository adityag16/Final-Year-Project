function [coord_region] = findRegionLW(array,start_pt)
%FOLLOWCONTOUR follows the contour of an array region of value 1 defining 
% by bounds.

% Constants
[ROWS,COLS,P] = getConstantsLW;

pts_to_check = [start_pt];
pts_visited = zeros(ROWS,COLS);
region = [];
while ~isempty(pts_to_check)
    N = size(pts_to_check,1); % number of points
    for i = 1:N
        pt = pts_to_check(i,:); % Current point
        % Extract subarray of neighbours of current point
        C_cross = [0 1 0; 1 0 1; 0 1 0];
        C_square = [1 1 1; 1 0 1; 1 1 1];
        M = zeros(size(array)); M(pt(1),pt(2))=1;
        n_cross = length(find(array(conv2(M,C_cross,'same')>0)<2));
        n_square = length(find(array(conv2(M,C_square,'same')>0)<2));

        
        % Decide whether to add current pixel to the region
        if (n_cross>1 && n_square>3) || isequal(pt,start_pt)
            region = [region;pt];
        end
        pts_visited(pt(1),pt(2))=1;
        
        % Iterate over neighbouring pixels
        % Cross-neighbours of current pixel
        [neighbours_x,neighbours_y] = find(conv2(M,C_cross,'same')==1);
        neighbours = [neighbours_x,neighbours_y];
        
        for j = 1:length(neighbours)
            % Add to list of points to check if not already there
            pt_n = [neighbours(j,1),neighbours(j,2)];
            if ~ismember(pt_n,region,'rows') && ~ismember(pt_n,pts_to_check,'rows') && ...
                pts_visited(pt_n(1),pt_n(2))==0
                if array(pt_n(1),pt_n(2))<2
                    pts_to_check = [pts_to_check;neighbours(j,:)];
                end
            end
        end
    end
    % Empty points to check
    pts_to_check(1:N,:) = [];
end

coord_region = sort(a2v(region,ROWS,COLS));
% Debug
% figclf(21);
% vect_out = zeros(1,ROWS*COLS);
% vect_out(coord_region) = 1;
% array_out = reshape(vect_out,ROWS,COLS);
% surf(array_out); view(2);