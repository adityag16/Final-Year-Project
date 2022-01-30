function out_smooth = smoothPixelsLW(out,smooth_factor,pixel_on)
%smoothPixelsLW filters the outputs of all active pixels

out_smooth = zeros(size(out));
for p = 1:npon
    out_smooth(:,pixel_on(p)) = ...
        smooth(out(:,pixel_on(p)),smooth_factor)-smooth(mean(out(:,coord_temp),2),smooth_factor);
    % Here add subtraction to closest temp pixel
    
end

end

