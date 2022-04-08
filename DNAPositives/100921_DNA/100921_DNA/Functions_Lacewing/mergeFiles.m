function cell_out = mergeFiles(filename_in,filename_out)
%mergeFiles merges two files containing cells

% Constants
N = length(filename_in);

% Initialisation
data = cells(N,1);

% Execution
for i = 1:N
    data{i} = importdata(filename_in{i});
    str = cellfun(@x,x.i);
end

end

