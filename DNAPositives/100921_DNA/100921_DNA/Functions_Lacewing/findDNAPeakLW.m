function [xm tm] = findDNAPeakLW(out,time,smooth_factor,range)
%findDNAPeakLW finds the max derivative of the sensor output xm

x_stable = range(1);
[~,~,np] = getConstantsLW;

P = size(range,1);
if nargin<4 % If no x0 specified, peak can be anywhere
    range_x = ones(P,1)*[1 length(time)];
else % If x0 specified, define an interval to find the peak
    range_x = [max([range(:,1),x_stable*ones(P,1)],[],2) min([range(:,2),length(time)*ones(P,1)],[],2)];
end

% Initialisation
Pon = size(out,1);
xm = zeros(Pon,1);

for p=1:Pon  
    % Derivative
    diff1_out = smooth(diff(out(p,:)),smooth_factor);
    % Uncomment to plot
%     figure(7); plot(time(2:end),diff1_out);
    
    % Find amplification peak xm
    [y,x,w] = findpeaks(diff1_out(range_x(p,1):range_x(p,2)-1),'Annotate','extents');
    if ~isempty(x)
        [~,idx] = max(w.*(y-min(y)));
        % Maximise a factor depending on the width of the peak and its
        % value. Y can be positive or negative so it is levelled up to the
        % positive values within the factor.
        xm(p) = x(idx)+range_x(p,1)-1;
    else
        [~,xm(p)] = max(diff1_out);
    end
end

tm = time(xm);

end
