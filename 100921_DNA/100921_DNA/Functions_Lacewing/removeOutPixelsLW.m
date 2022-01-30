function [pixel_on, pixel_off, npon,npoff] = removeOutPixelsLW(out,pixel_on,x_stable)
%removeOutPixels selects pixels that are in range for the whole experiment 

[~,~,P,coord_chem] = getConstantsLW;

if nargin<3
    x_stable = 1;
end
if nargin<2
    pixel_on = 1:P;
end

npon = length(pixel_on);
pixel_in_bounds = zeros(npon,1);
for p = 1:npon
    pixel_in_bounds(pixel_on(p)) = all(out(pixel_on(p),x_stable:end) >= 50 & out(pixel_on(p),x_stable:end) <= 950);
end
pixel_on = mintersect(find(pixel_in_bounds==1),coord_chem,pixel_on);
pixel_off = setdiff(1:P,pixel_on);
npon = length(pixel_on); % Number of on pixels
npoff = length(pixel_off); % Number of off pixels

end

