function array = coord_to_array(coord,X,Y, plot)
%Converts the coordinates to an array of dimension X x Y
% Also plots it if 1 is specified

if nargin<4
    plot = 0;
end

vect = zeros(X*Y,1);
vect(coord) = 1;
array = reshape(vect,X,Y);

if plot
    figure;surf(array); view(2);
end

end

