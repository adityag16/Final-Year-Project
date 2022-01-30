function [coord_wells,wells] = divideArrayLW(thr)
%divideArrayLW divides the array across regions determined by row thr, into
% length(thr)+1 region.

% Initialisation
[ROWS,COLS,np] = getConstantsLW;
wells = zeros(np,1);
coord_array = zeros(np,2);

% First generate all coordinates from the array
i=1;
for r = 1:ROWS
    for c = 1:COLS
        coord_array(i,:)=[r c];
        i = i+1;
    end
end

% Then convert to coordinates and check which region it belongs to
W = length(thr);
coord_vect = zeros(np,1);
for i = 1:np
    coord_vect(i) =a2v(coord_array(i,:),ROWS,COLS);
    wells(i) = W+1;
    for w = 1:W
        if coord_array(i,1) <= thr(w)
            wells(i) = w;
        end
    end
end

% Output coordinates
coord_wells=cell(W+1,1);
for w = 1:W+1
    coord_wells{w} = coord_vect(wells==w)';
end