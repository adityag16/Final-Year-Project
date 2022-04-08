function [h,array] = surfArrayLW(array_in,options,idx,n_fig)
%surfArrayLW plots the array

if nargin<2
    options = 'none';
end

[ROWS,COLS,np] = getConstantsLW;

switch options
    case 'none' % Just plot
        array = array_in;
    case 'reshape' % Just reshape to array
        array = reshape(array_in,ROWS,COLS);
    case 'fill' % Reshape and fill with idx
        vect = zeros(1,np);
        vect(idx) = array_in;
        array = reshape(vect,ROWS,COLS);
end

if nargin==4
    h = figure(n_fig);
end
surf(array); 
view(2);
colorbar;
axis equal;
xlim([1 COLS]); ylim([1 ROWS]);

% Scale  
vect = reshape(array,np,1);
vect = nonzeros(vect);
vect_noout = rmoutliers(vect);
% if std(vect_noout)>0
%     caxis([mean(vect_noout) - std(vect_noout) mean(vect_noout)+std(vect_noout)]);
% end
end

