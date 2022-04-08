function sigm_param = extractSigmParamLW(sigmoid_fit,sigmoid_scale,time,pH_init)
%extractSigmParamLW extract features of the sigmoidal amplification curve
% including Ct, efficiency, amplitude and asymmetry

% 4 figures for each relevant parameter
sig_5param = @(Fb, Fmax, inflection, slope, d, x) Fb + Fmax./(1 + exp(-(x - inflection)./slope)).^(1/d);
reverse = @(Fb, Fmax, inflection, slope, d, y) inflection - slope.*log( (Fmax./(y - Fb)).^d - 1 );
const = @(x, inflection, slope) (exp(inflection-x)./slope);
sigm5_deriv = @(Fmax, inflection, slope, d, x) (Fmax.*const(x, inflection, slope).*...
    (const(x, inflection, slope)+1).^(-(1+d)/d))./(d*slope);
% Parameters
Ct_thr = .1;

% Initialisation
sigm_param = zeros(1,4);
% Ct
Ct = time(end)*reverse(sigmoid_fit.Fb,sigmoid_fit.Fmax...
    ,sigmoid_fit.inflection,sigmoid_fit.slope,sigmoid_fit.d,Ct_thr)/60;
% Efficiency - slope
slope = max(sigm5_deriv(sigmoid_fit.Fmax...
    ,sigmoid_fit.inflection,sigmoid_fit.slope,sigmoid_fit.d,time));
% Amplitude in pH
ampl = sigmoid_fit.Fmax*(sigmoid_scale(2)-sigmoid_scale(1));
ampl = -log10((10^(-pH_init)+ampl)/10^(-pH_init));
% Asymmetry
asym = sigmoid_fit.d;   
% If parameters are real store them
if isreal(Ct) && isreal(slope) && isreal(ampl) && isreal(asym)
    sigm_param = [Ct slope ampl asym];
else
    [Ct slope ampl asym]
%                     figure(3); plot(time,sigmoid_fit(time))
    %sigm_param{i}{j}(p,:) = [0 0 0 0];
%                     sig_5param(sigmoid_fit.Fb,sigmoid_fit.Fmax...
%                         ,sigmoid_fit.inflection,sigmoid_fit.slope,sigmoid_fit.d,time)
end



end

