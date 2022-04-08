function [array_filter vect_filter] = movingAverageLW(vect, pixel_on)
%MOVINGAVERAGELW

	Pon= length(pixel_on);
    N = size(vect,2);
    
    % Import constants
    [ROWS,COLS,P,coord_chem,ROWS_temp,COLS_temp,npt,coord_temp,coord_peltier] = getConstantsLW;
    pixel_off = setdiff(1:P,pixel_on);
    
    % Convert vect into an array
    vect_all = zeros(P,N);
    vect_all(pixel_on,:) = vect;
        
    % Replace off pixels by NaN
    vect_all(pixel_off,:) = NaN;
    
    % Reshape vect to array
    array_all = reshape(vect_all,ROWS,COLS,N);
    array_filter = movmean(array_all,3,'omitnan');
    vect_filter = reshape(array_filter,P,N);
    vect_filter = vect_filter(pixel_on,:);

end