function [out_exp,out_compensated,rsquare] = compensateExp2LW(out,time,range,x_stable)
%compensateExp2LW fits an exponential to the start of the sensor output and
%   subtract it to the measurements

% Fitting options
alg='trust-region-reflective' ; % 'trust-region-reflective' or 'levenberg-marquardt'
displaySetting='off'; % 'iter', 'off', 'final', etc.
maxIter=1e3; % for debugging
maxFunEvals=1e3;
tolfun=1e-10;

% Initialisation
npon = size(out,2);
% out_exp = zeros(n,npon);
% out_compensated = zeros(n,npon);
% rsquare = zeros(npon,1);

optsN=optimoptions('lsqcurvefit','display',displaySetting,'Algorithm',alg);


range_wide = [min(range(:,1)),max(range(:,2))];
n = diff(range_wide)+1;
% [xx_mat,yy_mat] = initialiseArray([diff(range_wide)+1,npon]);
xx_mat=time(range_wide(1):range_wide(2))*ones(1,npon);
yy_mat = out(range_wide(1):range_wide(2),:);
range = range-range(1)+1;
range_wide = range_wide-range_wide(1)+1;

for p = 1:npon
    % Normalisation
    xx_mat(:,p) = (xx_mat(:,p) - xx_mat(range(p,1)))/...
        (xx_mat(range(p,2))-xx_mat(range(p,1)));
    yy_mat(:,p) = (yy_mat(:,p) - yy_mat(range(p,1),p))/...
        (yy_mat(range(p,2),p)-yy_mat(range(p,1),p));
    
%     figclf(11); plot(time(range_wide(p,1):range_wide(p,2)),yy_mat(:,p))
    
    % Correct start
    if range(p,1)>range_wide(1)
        diff_range = range(p,1)-range_wide(1)+1;
        yy_mat(x_stable:range(p,1),p) = yy_mat(range(p,1)+1:range(p,1)+diff_range+1,p); 
    end
    % Replace end of interpolating vector by a repeat of the relevant data
    % points
    if range(p,2)<range_wide(2)
        n_data_to_remove = range_wide(2) - range(p,2) + 1;
        n_data_to_repeat = diff(range(p,:));
        % Find how many filling iterations are necessary
        n_iter = ceil(n_data_to_remove/n_data_to_repeat);
        for i = 1:n_iter
            yy_mat(end-i*n_data_to_repeat+1:end-(i-1)*n_data_to_repeat,p) = ...
                yy_mat(range(p,2)-n_data_to_repeat+1:range(p,2),p);
            xx_mat(end-i*n_data_to_repeat+1:end-(i-1)*n_data_to_repeat,p) = ...
                xx_mat(range(p,2)-n_data_to_repeat+1:range(p,2),p);
        end 
    end
%     figclf(12); plot(time(range_wide(p,1):range_wide(p,2)),yy_mat(:,p))
end
c0 = ones(1,npon)'*[1000,1000,1];

% % Fitting test
% c0 = c0(1:5,:);
% xx_mat = xx_mat(:,1:5);
% yy_mat = yy_mat(:,1:5);

tic
[pfit]=lsqcurvefit(@fun_exp_drift,c0,xx_mat',yy_mat',[],[],optsN);
toc

% figclf(10);
% for p = 1:5
%     subplot(3,2,p); 
%     hold all;
%     plot(xx_mat(range(p,1):range(p,2),p),fun_exp_drift(pfit(p,:),xx_mat(range(p,1):range(p,2),p)));
%     plot(xx_mat(range(p,1):range(p,2),p),yy_mat(range(p,1):range(p,2),p));
% end

p

% Fitting
try
    tic
%         [f_opt,gof] = fit(x_opt(range(1):range(2)),y_opt(range(1):range(2)),exp2,options);
    f_opt = lsqcurvefit(exp2,[4e3, 3e4,1 ],x_opt(range(1):range(2)),y_opt(range(1):range(2)));
%         if abs(gof.rsquare)>0.2
%             range(2) = max(range(1)+10,round(range(2)*3/4));
%             [f_opt,gof] = fit(x_opt(range(1):range(2)),y_opt(range(1):range(2)),exp2,options);
%         end
    toc
    figure (10);
    plot(x_opt,y_opt,x_opt,exp2(f_opt,x_opt))
    out_exp(:,pixel_on(p)) = f_opt(x_opt)*(y(end)-y(x_stable))+y(x_stable);
    out_exp(1:x_stable-1,pixel_on(p)) = out_exp(x_stable,pixel_on(p));
    % Because sometimes weird problems of complex numbers
    out_compensated(:,pixel_on(p)) = (out(:,pixel_on(p)) - out_exp(:,pixel_on(p)));
    rsquare(pixel_on(p)) = gof.rsquare;
catch
    out_exp(:,pixel_on(p)) = zeros(n,1);
    out_compensated(:,pixel_on(p)) = out(:,pixel_on(p));
    rsquare(pixel_on(p)) = 0;
end



end

function [f]=fun_exp_drift(x,xdata)
    f = x(:,1).*(1-exp(-(xdata./x(:,2)).^x(:,3)));
end