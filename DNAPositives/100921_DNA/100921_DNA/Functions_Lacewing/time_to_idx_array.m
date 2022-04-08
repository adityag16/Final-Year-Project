function arr = idx_to_array(idx,)

arr = zeros(length(time), 1);

for i = 1:length(time)
    [~, arr(i)] = min(abs(time_vect-time(i)));

end

end

