function [xi,xf,out_exp_i,out_exp_f] = findDNABoundsLW(out,time,x_stable,xm,smooth_factor)
%findDNABounds find the bounds of the amplification for the sensor output

% Initialisation
N = size(out,2);
Pon = size(out,2);

% Linear interpolation around xm
pol = zeros(N,2);
y_pol = zeros(size(out));
for p = 1:Pon
    pol(p,:) = polyfit(time(xm(p)-5:xm(p)+5),out(p,xm(p)-5:xm(p)+5),1);
    y_pol(:,p) = polyval(pol(p,:),time');
end

%  Find deviation points
diff1_out = zeros(N-1,Pon);
[xi,xf] = initialiseArray([Pon,2]);
t_i = time_to_idx(time,time(x_stable) + [1 2]*60);
t_f = time_to_idx(time,[time(end)-3*60 time(end)-2*60]);
for p= 1:Pon
    diff1_out(p,:) = smooth(diff(smooth(out(p,:),smooth_factor)),3*smooth_factor);
    
    % Debug - plot
    if p>1000
        p
        figclf(11);
        h(1)=subplot(2,1,1); hold all;
        plot(time,out(p,:));
        plotValue(xm(p),time,out(p,:),'X');
        h(2)=subplot(2,1,2); hold all;
        plot(time(2:end),diff1_out(p,:));
        plotValue(xm(p),time,diff1_out(p,:),'X');
        linkaxes(h,'x')
        p;
    end
    
    diff1_xi = mean(diff1_out(p,t_i(1):t_i(2)));
    % Find dip after xm
    [dip_y,dip_x] = findpeaks(-diff1_out(p,:));
    idx_dip = intersect(find(dip_x>xm(p)),find(-dip_y<diff1_out(p,xm(p))));
    if ~isempty(idx_dip)
        diff1_xf = diff1_out(p,dip_x(idx_dip(1)));
    else
        diff1_xf = min(diff1_out(p,xm(p):end-2*60));
        if isempty(diff1_xf)
            diff1_xf = min(diff1_out(p,xm(p):end));
        end
    end
    diff1_xf_2 = min(diff1_out(p,xm:end-10));
    
    thr_i = .9;
    thr_f = .3;
    dist_x1 = diff1_out(p,xm(p))-(1-thr_i)*abs(diff1_out(p,xm(p))-diff1_xi); % Distance threshold
    dist_x2 = diff1_out(p,xm(p))-(1-thr_f)*abs(diff1_out(p,xm(p))-diff1_xf);
    dist_x2_2 = diff1_out(p,xm(p))-(1-thr_f)*abs(diff1_out(p,xm(p))-diff1_xf_2);
    x1 = find(diff1_out(p,x_stable:xm(p))-dist_x1<0);
    x2 = find(diff1_out(p,xm(p):end)-dist_x2<0);
    x2_2 = find(diff1_out(p,xm(p):end)-dist_x2_2<0);
    if ~isempty(x1)
        x1 = x_stable+x1(end)-1;
    else
        x1=xm(p)-2;
    end
    if ~isempty(x2)
        x2 = x2(1);
    else
        x2= 2;
    end
    x2 = min(x2 + xm(p)-1,length(time)); % Leveling up
    if ~isempty(x2_2)
        x2_2 = x2_2(1);
    else
        x2_2= 2;
    end
    x2_2 = min(x2_2 + xm(p)-1,length(time)); % Leveling up
    
    % Fill matrices
    xi(p,1) = x1; xf(p,1)= x2; xf(p,2) = x2_2;
end


% Adding a min on length time to make sure there is no overflowing

% [~,x1] = min(abs(diff1_out(x_stable:xm)-dist_x1));
% [~,x2] = min(abs(diff1_out(xm:end)-dist_x2));

% Exponential fitting
% Two types of initial points and lowest error between both

% Initial fitting close to x1
range_i(:,:,1) = [x_stable*ones(Pon,1),xi(:,1)];
disp('Interpolation 1')
[out_exp_i(:,:,1),res_i(:,1),pfit_i(:,:,1),convex_i(:,1)] = compensateExp2ParallelLW(out,time,range_i(:,:,1));
% Check if another method is better
for p = 1:Pon
    [~,idx] = min(abs(out_exp_i(p,x_stable:xm)-y_pol(x_stable:xm,p)));
    xi(p,2) = idx(1) + x_stable - 1;
end
range_i(:,:,2) = [x_stable*ones(Pon,1),xi(:,2)];
disp('Interpolation 2')
[out_exp_i(:,:,2),res_i(:,2),pfit_i(:,:,2),convex_i(:,2)] = compensateExp2ParallelLW(out,time,range_i(:,:,2));

% Final fitting close to x2
range_f(:,:,1) = [xf(:,1),(length(time)-10)*ones(Pon,1)];
range_f(:,:,2) = [xf(:,2),min(xf(:,2)+200,length(time))];

disp('Interpolation 3')
convex_in = sign(mean(convex_i(:,1)));

[out_exp_f,res_f,pfit_f,convex_f,range_f] = compensateExp2ParallelLW(out,time,range_f,pfit_i(:,:,2),convex_in);

% Plot all
for p= 1:Pon
    if p>1000
        
        figclf(12);
        subplot(2,2,1); hold all;
        plot(time,out(:,p),'--');
        plot(time,out_exp_i(p,:,1));
        plot(time,out_exp_f(p,:,1));
    %     title(['Residual is ' num2str(res_i_1(p))])
        ylim([min(out(p,:)) max(out(p,:))])
        title(['Convexity = ' num2str(convex_i(p,1))])

        subplot(2,2,2); hold all;
        plot(time,out(p,:),'--');
        plot(time,out_exp_i(p,:,1));
        plot(time,out_exp_f(p,:,2));
        ylim([min(out(p,:)) max(out(p,:))])
        title(['Convexity = ' num2str(convex_i(p,1))])

        subplot(2,2,3); hold all;
        plot(time,out(p,:),'--');
        plot(time,out_exp_i(p,:,2));
        plot(time,out_exp_f(p,:,1));
        ylim([min(out(p,:)) max(out(p,:))])
        title(['Convexity = ' num2str(convex_i(p,2))])

        subplot(2,2,4); hold all;
        plot(time,out(p,:),'--');
        plot(time,out_exp_i(p,:,2));
        plot(time,out_exp_f(p,:,2));
    %     title(['Residual is ' num2str(res_i_2(p))])
        ylim([min(out(:,p)) max(out(p,:))]) 
        title(['Convexity = ' num2str(convex_i(p,2))])
    end
end

end

