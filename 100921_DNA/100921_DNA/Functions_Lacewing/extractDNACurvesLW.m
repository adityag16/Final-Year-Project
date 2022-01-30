function [x0,xi,xf,ampl_curve_V,exp_comp,ampl_curve_pH,ampl_curve_lin,ampl_curve_sigm,ampl_curve_norm,time_sigm,outcome,sigm_param,range] ...
    = extractDNACurvesLW(out,time,smooth_factor,range,pH_init,wb)
%extractDNACurves converts a single or a set of raw output voltage curves 
% to a single or a set of amplification curves.

x_stable = range(1);
Np = range(2);

% Constants
N = size(out,2);
P = size(out,1);
V = 4; % Number of variable parameters for sigmoidal fitting

% Initialisation
[out_exp_i,out_exp_f,ampl_curve_V,ampl_curve_pH,ampl_curve_lin,exp_comp,ampl_curve_sigm] = initialiseArray([P,N]);
[xi,xf,error_fit,outcome] = initialiseArray([P,1]);
sigm_param= zeros(V,P);

% Find preliminary maximum of derivative
[x0,t0] = findDNAPeakLW(out,time,smooth_factor,ones(P,1)*[x_stable Np]);

Np = min(Np,time_to_idx(time,mean(t0)+1200));

% Subtract initial fitting to get more accurate x0
% Range is decided from a minute after stable time to half of x0
range = ones(P,1)*time_to_idx(time,time(x_stable) + [60 mean(time(x0))/2]);
[out_exp_init,res_f] = compensateExp2ParallelLW(out,time,range);
out_comp = out - out_exp_init;

% % Debug - plot
% for p = 1001:P
%     figclf(15);
%     subplot(2,2,1); hold all;
%     plot(time,out(:,p));
%     plot(time,out_exp_init(:,p));
%     subplot(2,2,2);
%     plot(time,out_comp(:,p))
%     subplot(2,2,3); hold all;
%     plot(time(2:end),diff(out(:,p)));
%     plot(time(2:end),diff(out_exp_init(:,p)));
%     subplot(2,2,4);
%     plot(time(2:end),diff(out_comp(:,p)))
% end

range = [max(x_stable,time_to_idx(time,t0-12*60)') Np*ones(P,1)]; % Range for interpolation
% range = [x_stable Np]; % Range for interpolation

% Find accurate maximum of derivative
xm = findDNAPeakLW(out_comp,time,smooth_factor,range);

% Find amplification region for average
[out_exp_i,out_exp_f,xi,xf] = findDNABounds2LW(out,time,range,x_stable,xm);

% Drift compensation
[ampl_curve_V,exp_comp] = compensateDriftLW(time,out,xi,xf,out_exp_i,out_exp_f,range(1));

% Linearise amplification curves
[ampl_curve_lin,ampl_curve_pH] = linAmplCurveLW(ampl_curve_V,range(1),pH_init,smooth_factor);

% Sigmoidal fitting - Disables
[ampl_curve_sigm,sigmoid_fit,time_sigm,sigmoid_gof,sigmoid_scale] = fitSigmoidLW(ampl_curve_lin,time,range);
% ampl_curve_sigm = ampl_curve_lin;
% sigmoid_fit = 0;
% time_sigm = 0;
% sigmoid_gof = 0;
% sigmoid_scale = 0;
% sigm_param = 0;

% Normalisation
ampl_curve_norm = (ampl_curve_sigm-ampl_curve_sigm(:,1))./(ampl_curve_sigm(:,end)-ampl_curve_sigm(:,1));
% ampl_curve_norm(:,pixel_off) = [];


% error_fit = sigmoid_gof.rsquare;
% 
% % Sigmoidal feature extraction
% if error_fit(i)>.95 && (ampl_curve_lin(end,i)>ampl_curve_lin(1,i))
%     sigm_param(:,i) = extractSigmParamLW(sigmoid_fit,sigmoid_scale,time,pH_init);
%     outcome(i) = 1;
% else
%     outcome(i) = 2; % might be negative
% end
% 
% % Check negative
% if outcome(i) == 2
%     % Check if negative by fitting
%     [exp_comp(:,i),ampl_curve_V(:,i),rsquare_neg] = compensateExp2LW(out(:,i),time,[time(x_stable),time(end)],x_stable);
%     ampl_curve_V(1:x_stable,i) = 0;
%     if rsquare_neg<.95 % unstable output
%         outcome(i) = 0;
%     else
%         [ampl_curve_lin(:,i),ampl_curve_pH(:,i)] = linAmplCurveLW(ampl_curve_V(:,i),x_stable,pH_init,smooth_factor);
%         ampl_curve_sigm(:,i) = ampl_curve_lin(:,i);
%     end
% %         ampl_curve_V(:,i) = 500*ones(size(out(:,i)));
% %         ampl_curve_pH(:,i) = pH_init*ones(size(out(:,i)));
% %         ampl_curve_lin(:,i) = 10^(-pH_init)*ones(size(out(:,i)));
% %         ampl_curve_sigm(:,i) = ampl_curve_lin(:,i);
end