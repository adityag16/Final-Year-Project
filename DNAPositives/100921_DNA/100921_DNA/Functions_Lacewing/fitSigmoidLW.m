function [ampl_curve_sigm,sigmoid,time_sigm,sigmoid_gof,sigmoid_scale] = fitSigmoidLW(ampl_curve_lin,time,range)

% Constants
P = size(ampl_curve_lin,1);
N = size(ampl_curve_lin,2);

% Fitting options
options = fitoptions('Method', 'NonlinearLeastSquares',...
                     'Lower', [-Inf, -Inf, -Inf, -Inf, 0],...
					 'Upper', [Inf, Inf, Inf, Inf, 2],...
                     'Robust', 'Off'); %'BiSquare');

% Initialisation
time_sigm = time(1):1:time(end);
ampl_curve_sigm = zeros(P,length(time_sigm));
sigmoid_gof = cell(P,1);
sigmoid_scale = zeros(P,2);
for p = 1:P
    % Initialisation
    x = time;
    y = ampl_curve_lin(p,:);
    sig_4param = @(Fb, Fmax, inflection, slope, d, x) Fb + Fmax./(1 + exp(-(x - inflection)./slope)).^(1/d);
    options.StartPoint = [0, 1, 0.5, 0.04, 1];
    
    % Vector normalisation
    x_opt = x/x(end); 
    if y(end) ~= y(1)
        y_opt = (y - y(1))/(y(end) - y(1));
    else
        y_opt = y;
    end
    [sigmoid,gof] = fit(x_opt(range(p,1):range(p,2))',y_opt(range(p,1):range(p,2))',sig_4param,options);            
    sigmoid_gof{p} = gof;
    ampl_curve_sigm(p,:) = sigmoid(time_sigm/time_sigm(end))*(y(end)-y(1));
    ampl_curve_sigm(p,:) = ampl_curve_sigm(p,:) - ampl_curve_sigm(p,1)+y(1);
    sigmoid_scale(p,:) = [y(1) y(end)]; % Scaling factor
end

