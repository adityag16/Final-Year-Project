function [x1,x2] = findDNABounds(out,time,smooth_factor,x_stable,xm,pixel_on)
%findDNABounds find the bounds of the amplification for the sensor output

[~,~,np] = getConstantsLW;
if nargin<6
    if size(out,2)==1 % Only one curve
        pixel_on = 1;
    else % For all pixels
        pixel_on = 1:np;
    end
end

% Initialisation
npon = length(pixel_on);
x1 = zeros(npon,1);
x2 = zeros(npon,1);
for p = 1:npon
    if xm(pixel_on(p))~=0
        % Derivative
        diff1_out = smooth(diff(out(:,pixel_on(p))),smooth_factor);
        
        % Linear aproximation of xm
        pol_m = polyfit(time(xm(pixel_on(p))-10:xm(pixel_on(p))+10),...
            out(xm(pixel_on(p))-10:xm(pixel_on(p))+10,pixel_on(p))',1);
        y_pol_m(:,pixel_on(p)) = polyval(pol_m,time);
        
        %  Find deviation points
        x1_idx = find(abs(out(:,pixel_on(p))-y_pol_m(:,pixel_on(p)))>.5);
        idx = find(x1_idx<xm(pixel_on(p)));
        if isempty(idx) || x1_idx(idx(end)) < x_stable + 3
            % Bad pixel
            x1(pixel_on(p)) = 5;
            xm(pixel_on(p)) = 0;
            pixel_on_dna(pixel_on_dna==pixel_on(p)) = []; % Remove from dna pixels
        else
            x1(pixel_on(p)) = x1_idx(idx(end));
        end
        x2_idx = find(abs(out(:,pixel_on(p))-y_pol_m(:,pixel_on(p)))>.5);
        idx = find(x2_idx>xm(pixel_on(p)));
        if ~isempty(idx)
            x2(pixel_on(p)) = x2_idx(idx(1));
        else
            x2(pixel_on(p)) = n-4;
            xm(pixel_on(p)) = 0;
            pixel_on_dna(pixel_on_dna==pixel_on(p)) = []; % Remove from dna pixels
        end
        
        figclf(4); 
        h(1)=subplot(2,1,1); hold all; plot(time,out(:,pixel_on(p)));
        plot(time(xm(pixel_on(p)))*ones(1,100),linspace(min(out(:,pixel_on(p))),max(out(:,pixel_on(p))),100));
        plot(time(x1(pixel_on(p)))*ones(1,100),linspace(min(out(:,pixel_on(p))),max(out(:,pixel_on(p))),100));
        plot(time(x2(pixel_on(p)))*ones(1,100),linspace(min(out(:,pixel_on(p))),max(out(:,pixel_on(p))),100));
        h(2)=subplot(2,1,2); hold all; plot(time(1:end-1),diff1_out);
        plot(time(xm(pixel_on(p)))*ones(1,100),linspace(min(diff1_out),max(diff1_out),100));
        plot(time(x1(pixel_on(p)))*ones(1,100),linspace(min(diff1_out),max(diff1_out),100));
        plot(time(x2(pixel_on(p)))*ones(1,100),linspace(min(diff1_out),max(diff1_out),100));
        linkaxes(h,'x');
        %waitforbuttonpress;
    end

end

% Approximation of x1
t = 3:floor(time(xm)); % Interval between 3 mins and the amplification point
for i = 1:length(t)
    [out_exp,out_comp_1,rsquare(i)] = compensateExp2LW(out_smooth_mean,time,[time(end)-t(i)*60,time(end)]);
    subplot(4,4,i);
    plot(time/60,out_smooth_mean,time,out_exp);
    title(['t = ' num2str(t(i)) ' s'])
end