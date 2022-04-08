function idx_v = a2v(idx_a,ROWS,COLS)
%v2a converts vector coordinates to array coordinates

N = size(idx_a,1); % number of coordinates to convert
idx_v = zeros(N,1);
for i = 1:N
    idx_v(i) = (idx_a(i,2)-1)*ROWS+idx_a(i,1);
end

