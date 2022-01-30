function idx = time_to_idx(time_vect, time)

for i = 1:length(time)
    [~, idx(i)] = min(abs(time_vect-time(i)));
end

end

