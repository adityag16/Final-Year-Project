function [ampl_curve_lin,ampl_curve_pH] = linAmplCurveLW(ampl_curve_V,x_stable, pH_init,smooth_factor)
%linAmplCurveLW linerarises the chip output to produce an amplification
%curve following the concentration of protons in solution

% Constants
P = size(ampl_curve_V,1);
N = size(ampl_curve_V,2);
sens = 13; % mV/pH

% Initialisation
[ampl_curve_pH,ampl_curve_lin] = initialiseArray([P,N]);

for p = 1:P
    % Generate output amplification curve
    ampl_curve_pH(p,:) = smooth(ampl_curve_V(p,:)-ampl_curve_V(p,x_stable),smooth_factor);
    ampl_curve_pH(p,:) = -ampl_curve_pH(p,:)/sens+pH_init;
    ampl_curve_lin(p,:) = smooth(10.^(-ampl_curve_pH(p,:)),smooth_factor);
end

end

