function [xi,xf,out_exp_i,out_exp_f] = findDNABounds(out,time,x_stable,xm)
%findDNABounds find the bounds of the amplification for the sensor output

% Approximation of xi
% Initialisation
t = 3*60:60:floor(time(xm)); % Interval between 3 mins and the amplification point
rsquare_i = zeros(length(t),1);
out_exp_i = cell(length(t),1);
figclf(2);
for i = 1:length(t)
    [out_exp_i{i},~,rsquare_i(i)] = compensateExp2LW(out,time,[time(x_stable),t(i)]);
    subplot(4,4,i);
    plot(time/60,out,time/60,out_exp_init);
    title(['t = ' num2str(t(i)) ' s'])
end
[~,idx_i] = max(rsquare_i);
out_exp_i = out_exp_i{idx_i};

% x1 found based on deviation
thr = .1*abs(out(xm)-out_exp_i(xm));
idx = find(abs(out-out_exp_i)<thr);
xi = intersect(idx,1:xm);
xi=xi(end);

% [out_exp_init,out_comp_init,rsquare_init(i)] = compensateExp2LW(out,time,[time(end)-t(i)*60,time(end)]);

% Approximation of xf
% Initialisation
t = ceil(time(xm)):60:time(end)-3*60; % Interval the amplification point and 3 mins before the end
rsquare_f = zeros(length(t),1);
out_exp_f = cell(length(t),1);
figclf(2);
for i = 1:length(t)
    [out_exp_f{i},~,rsquare_f(i)] = compensateExp2LW(out,time,[t(i) time(end)]);
    subplot(5,5,i);
    plot(time/60,out,time/60,out_exp_f{i});
    title(['t = ' num2str(t(i)) ' s'])
end
figure(3); plot(t,rsquare_f);

[~,idx_f] = max(rsquare_f);
out_exp_f = out_exp_f{idx_f};

% x1 found based on deviation
thr = .1*abs(out(xm)-out_exp_f(xm));
idx = find(abs(out-out_exp_f)<thr);
xf = intersect(idx,xm:length(time));
xf=xf(1);

end

