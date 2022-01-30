function [ROWS,COLS,np,coord_chem,ROWS_temp,COLS_temp,npt,coord_temp,coord_peltier] = getConstantsLW
% Loads constants for TITANICKS array.

ROWS = 78;
COLS = 56;
np = ROWS*COLS;

ROWS_temp = ceil((ROWS-1)/3);
COLS_temp = ceil((COLS-1)/3);
npt = ROWS_temp*COLS_temp;
row_temp = [2:3:ROWS];
col_temp = [2:3:COLS];
coord_temp = zeros(1,npt);
i = 1;
for c = 1:COLS_temp
    for r = 1:ROWS_temp 
        coord_temp(i) = (col_temp(c)-1) * ROWS + row_temp(r);
        i = i+1;
    end
end
coord_chem = setdiff(1:np,coord_temp);

row_peltier = 29:3:50;
col_peltier = 32:3:56;
nptp = length(row_peltier)*length(col_peltier);
coord_peltier = zeros(1,nptp);
i=1;
for c = 1:length(col_peltier)
    for r = 1:length(row_peltier) 
        coord_peltier(i) = (col_peltier(c)-1) * ROWS + row_peltier(r);
        i = i+1;
    end
end

end

