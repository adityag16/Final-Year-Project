function [pixel_on,diff_Vref] = compensateVrefLW(out,dV)
%compensateVref selects all pixels responding to a Vref sweep of 100 mV
% Time stamps are correct for Lacewing.

if nargin<2
    dV = 20; % Expected attenuation of 0.7
end

diff_Vref = abs(mean([out(:,7)-out(:,6), out(:,9)-out(:,10)],2));
pixel_on = find(diff_Vref>dV);

end

