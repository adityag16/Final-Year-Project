function [] = rangeCorrection()
%RANGECORRECTION performs correction to 

function [out_exp,res,pfit,convex_out] = compensateExp2ParallelLW(out,time,range_dim,options)
%compensateExp2LW fits an exponential to the start of the sensor output and
%   subtract it to the measurements
% - out (S x F) -> (N x F) defined by range
% - time (S x 1) -> (N x 1) defined by range
% - range_dim (F x 2 x H)

% Fitting parameters
S = length(time);
F = size(out,2); % Number of functions
P = 3; % Number of parameters
M = 1;% Number of multiples
H = size(range_dim,3); % Number of dimensions

% Options
if nargin<4
    options = '';
    convex_in = 0;
else
    [convex_in] = cellparam(options,{'convex'});
end

% Initialisation
out_exp = zeros(S,F,H);
res = zeros(F,H);
pfit = zeros(F,P,H); 
convex_out = zeros(F,H);

for h = 1:H
    % Range definition
    range = range_dim(:,:,h);
    range_wide = [min(range(:,1)),max(range(:,2))];

    % Fitting parameters
    N = diff(range_wide)+1;

    % Parameters
    alg='levenberg-marquardt' ; % 'trust-region-reflective' or 'levenberg-marquardt'
    displaySetting='off'; % 'iter', 'off', 'final', etc.
    tol=5e-3;
    optsN=optimoptions('lsqcurvefit','display',displaySetting,'Algorithm',alg,'FunctionTolerance',tol);
    warning('on')

    % Data generation
    xx_mat=time(range_wide(1):range_wide(2))*ones(1,F);
    yy_mat = out(range_wide(1):range_wide(2),:);

    % Range offset
    range_corr = range-range_wide(1)+1;
    range_wide_corr = range_wide-range_wide(1)+1;

    % Normalisation
    % x
    idx = range_corr+[0:N:N*(F-1)]'*ones(1,2);
    xx_mat_norm_coeff = xx_mat(idx);
    xx_mat_norm = (xx_mat - ones(N,1)*xx_mat_norm_coeff(:,1)')./...
        (ones(N,1)*xx_mat_norm_coeff(:,2)'-ones(N,1)*xx_mat_norm_coeff(:,1)');
    x_norm =(time*ones(1,F) - ones(S,1)*xx_mat_norm_coeff(:,1)')./...
        (ones(S,1)*xx_mat_norm_coeff(:,2)'-ones(S,1)*xx_mat_norm_coeff(:,1)');
    % y
    yy_mat_norm_coeff = yy_mat(idx);
    yy_mat_norm = (yy_mat - ones(N,1)*yy_mat_norm_coeff(:,1)')./...
        (ones(N,1)*yy_mat_norm_coeff(:,2)'-ones(N,1)*yy_mat_norm_coeff(:,1)');

    % Range correction
    xx_mat_norm_corr = xx_mat_norm;
    yy_mat_norm_corr = yy_mat_norm;
    for f = 1:F
        % Correct start
        if range(f,1)>range_wide(1)
            n_data_to_remove = range_corr(f,1)-1;
            n_data_to_repeat = diff(range_corr(f,:))+1;
            % Find how many filling iterations are necessary
            n_iter = ceil(n_data_to_remove/n_data_to_repeat);
            for i = 1:n_iter
                range_remove = max(n_data_to_remove-i*n_data_to_repeat+1,1):n_data_to_remove-(i-1)*n_data_to_repeat;
                range_repeat = range_corr(f,1):range_corr(f,1)+length(range_remove)-1;
                yy_mat_norm_corr(range_remove,f) = yy_mat_norm(range_repeat,f);
                xx_mat_norm_corr(range_remove,f) = xx_mat_norm(range_repeat,f);
            end 
        end
        % Replace end of interpolating vector by a repeat of the relevant data
        % points
        if range(f,2)<range_wide(2)
            n_data_to_remove = range_wide_corr(2) - range_corr(f,2) + 1;
            n_data_to_repeat = diff(range_corr(f,:))+1;
            % Find how many filling iterations are necessary
            n_iter = ceil(n_data_to_remove/n_data_to_repeat);
            for i = 1:n_iter
                range_remove = i*n_data_to_repeat+1:min((i+1)*n_data_to_repeat,range_wide_corr(2));
                range_repeat = range_corr(f,2)-length(range_remove)+1:range_corr(f,2);
                yy_mat_norm_corr(range_remove,f) = yy_mat_norm(range_repeat,f);
                xx_mat_norm_corr(range_remove,f) = xx_mat_norm(range_repeat,f);
            end 
        end
    end
    p0 = ones(1,F)'*[150,650,.8];
%     p0 = ones(1,F)'*[-1,-1,1,0.5];

    tic;
    I = ceil(F/M);
    for i = 1:I
        i
        range_m = (i-1)*M+1:min(i*M,F);
        % If concavity option on, converge until the right concavity is
        % reached
        while (convex_out(range_m,h) == 0) || (convex_in ~=0 && convex_out(range_m,h) ~= convex_in && range(range_m(:),1) <= range(range_m(:),2) - 10)
            % Cycle through functions where concavity is wrong
            if (convex_out(range_m,h) ~=0 && convex_out(range_m,h) ~= convex_in)
                if i>1000
                    i
                    figclf(14); hold all;
                    plot(time,out(:,range_m));
                    plot(time,out_exp(:,range_m,h));
                    end
                idx = find(convex_out(range_m,h) ~= convex_in);
                for j = 1:length(idx)
                    range(range_m(idx(j)),1) = range(range_m(idx(j)),1) + 40;
%                     title(['Pixel ' num2str(i) ' range ' num2str(range(range_m(idx(j)),1))])
                end
            end
            
            % Fitting
            try
                [pfit(range_m,:,h),res(range_m,h)]=lsqcurvefit(@fun_exp_drift,p0(range_m,:),xx_mat_norm_corr(1:10:N,range_m)',yy_mat_norm_corr(1:10:N,range_m)',[],[],optsN);
            catch
                pfit(range_m,:,h) = p0(range_m,:);
                res(range_m,h) = 1e3*ones(M,1);
                disp('Fitting error')
            end
            
            % Function recovery
            for j = 1:length(range_m)
                range_itv = range(range_m(j),1):range(range_m(j),2);
                range_corr_itv = range_corr(range_m(j),1):range_corr(range_m(j),2);
                out_exp(:,range_m(j),h) = fun_exp_drift(pfit(range_m(j),:,h),x_norm(:,range_m(j))')'...
                    .*(yy_mat_norm_coeff(range_m(j),2)'-yy_mat_norm_coeff(range_m(j),1)')+...
                    yy_mat_norm_coeff(range_m(j),1)';
                out_diff2 = diff(diff(abs(out_exp(:,range_m(j),h))));
                convex_out(range_m(j),h) = sign(mean(out_diff2(min(range_itv,S)-2,:)));
            end
        end
    %     % Plot results
    %     figclf(10); 
    %     for j = 1:M
    %         subplot(M,1,j); hold all;
    %         plot(xx_mat_norm(:,range_m(j)),yy_mat_norm(:,range_m(j)),'-');
    %         plot(xx_mat_norm(:,range_m(j)),fun_exp_drift(pfit(range_m(j),:),xx_mat_norm(:,range_m(j))')','--');
    %     end
    %     i
    end
    toc;
end

end

function [f]=fun_exp_drift(x,xdata)
    f = x(:,1).*(1-exp(-(xdata./x(:,2)).^x(:,3)));
%     f = x(:,1).*exp(x(:,2).*xdata) + x(:,3).*exp(x(:,4).*xdata);
end
end

