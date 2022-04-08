function array_filter = movingAverage(array)
% Computes the 3-neighbour moving average of the array.

    dim = size(array);

    array = [array(1,:);array;array(end,:)];
    array = [array(:,1) array array(:,end)];
    array_filter = filter2(fspecial('average',3),array);
    array_filter = array_filter(2:dim(1)+1,2:dim(2)+1);
end