function array = pixel_on_to_pixel(coord,X,Y, pixel_on, plot)
%Converts the coordinates to an array of dimension X x Y
% Also plots it if 1 is specified

if nargin<5
    plot = 0;
end

vect = zeros(X*Y,1);
vect(pixel_on) = coord;
array = reshape(vect,X,Y);

if plot
    figure;surf(array); view(2);
    colorbar;
    axis equal;
    xlim([1 Y]); ylim([1 X]);
end

end

