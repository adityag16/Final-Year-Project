function [f]=fun_exp_drift(x,xdata)
    f = x(:,1).*(1-exp(-(x_data./x(:,2)).^x(:,3)));
end