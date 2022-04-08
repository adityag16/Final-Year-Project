function [ampl_curve_V,exp_comp] = compensateDriftLW(time,out,xi,xf,out_exp_i,out_exp_f,x_stable)
%compensateDriftLW removes drift based on previous exponential
%   interpolation

% Constants
P = size(out,1);
N = size(out,2);

% Initialisation
exp_comp = zeros(P,N);
ampl_curve_V = zeros(P,N);
for p = 1:P
    % Compensation curve
    exp_comp(p,:) = out_exp_i(p,:);
    % For the middle interval, calculate start and end slope
    if xi(p) > 7 && xf(p) < N-6
        slope_start = (out_exp_i(p,xi(p)-5)-out_exp_i(p,xi(p)-6))/(time(xi(p)-5)-time(xi(p)-6));
        slope_end = (out_exp_f(p,xf(p)+6)-out_exp_f(p,xf(p)+5))/(time(xf(p)+6)-time(xf(p)+5));
        drift_comp = (slope_end-slope_start)/(time(xf(p))-time(xi(p)))...
        *(time(xi(p):xf(p))-time(xi(p)))+slope_start;
        for j = xi(p)+1:xf(p)
            exp_comp(p,j) = exp_comp(p,j-1) + ...
                drift_comp(j-xi(p))*(time(j)-time(j-1));
        end
        exp_comp(p,xf(p):N) = out_exp_f(p,xf(p):N)-out_exp_f(p,xf(p)) + exp_comp(p,xf(p));
    end

    % Compensated voltage output
    ampl_curve_V(p,:) = out(p,:) - exp_comp(p,:) + out(p,x_stable) ;
%     ampl_curve_V(p,:) = ampl_curve_V(p,:) - ampl_curve_V(p,x_stable) + 500;
%     ampl_curve_V(p,1:x_stable-1) = 500;
    ampl_curve_V(p,1:x_stable-1) = ampl_curve_V(p,x_stable);
end

