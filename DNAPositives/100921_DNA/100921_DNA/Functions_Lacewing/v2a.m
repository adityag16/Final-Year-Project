function idx_a = v2a(idx_v,ROWS,COLS)
%v2a converts vector coordinates to array coordinates

N = length(idx_v); % number of coordinates to convert
idx_a = zeros(N,2);
for i = 1:N
    idx_a(i,1) = mod(idx_v(i)-1,ROWS)+1;
    idx_a(i,2) = ceil(idx_v(i)/ROWS);
end
    
end

