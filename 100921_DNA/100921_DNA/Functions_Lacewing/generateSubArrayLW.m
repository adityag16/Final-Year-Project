function coord = generateSubArrayLW(size)
%generateSubArrayLW returns the coordinates for a sub-array of a given size

    % Import LW constants
    [ROWS,COLS] = getConstantsLW;

    % Generate row and col
    d_row = floor(ROWS/size(1));
    d_col = floor(COLS/size(2));
    row = floor(d_row/2):d_row:floor(ROWS-d_row/2);
    col = floor(d_col/2):d_col:floor(COLS-d_col/2);
    % Convert to array coordinates
    coord = zeros(prod(size),1);
    for i = 1:prod(size)
        j = mod(i-1,size(1))+1;
        k = ceil((i/size(1)));
        coord(i) = (col(k)-1)*ROWS+row(j);
    end
    npon = length(coord);
end

