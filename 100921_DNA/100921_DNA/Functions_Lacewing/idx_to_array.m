function array = idx_to_array(vect,idx,ROWS,COLS)

array_vect=zeros(1,ROWS*COLS);
array_vect(idx) = vect;
array = reshape(array_vect,ROWS,COLS);

end

