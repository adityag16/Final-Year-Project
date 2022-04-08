function [out_exp,res,pfit,convex_out,range] = compensateExp2ParallelLW(out,time,range_dim,p0,dx,tol,convex_in)
%compensateExp2LW fits an exponential to the start of the sensor output and
%   subtract it to the measurements
% - out (F x S) -> (F x N) defined by range
% - time (1 x S) -> (1 x N) defined by range
% - range_dim (F x 2 x H)

% Fitting parameters
N = length(time);
F = size(out,1); % Number of functions
V = 3; % Number of variables/parameters
H = size(range_dim,3); % Number of dimensions

% Options
% Initial parameters for fitting 
if nargin<4 || isempty(p0)
    p0 = repmat([150,650,.8],[F,1,H]);
%     p0 = ones(1,F)'*[-1,-1,1,0.5];
end

% Sample increment
if nargin<5
    dx=10;
end

% Tolerance
if nargin<6
    tol=5e-3;
end

% Convexity restriction
if nargin<7
    convex_in = 0;
end

% Initialisation
out_exp = zeros(F,N,H);
res = zeros(F,H);
pfit = zeros(F,V,H); 
convex_out = zeros(F,H);

for h = 1:H
    % Range definition
    range = range_dim(:,:,h);

    % Parameters
    alg='levenberg-marquardt' ; % 'trust-region-reflective' or 'levenberg-marquardt'
    displaySetting='off'; % 'iter', 'off', 'final', etc.
    optsN=optimoptions('lsqcurvefit','display',displaySetting,'Algorithm',alg,'FunctionTolerance',tol);
    warning('on')

    % Data generation
    xx_mat=ones(F,1)*time;
    yy_mat = out;

    xx_mat_norm = zeros(F,N);
    yy_mat_norm = zeros(F,N);

%     tic;
    for f = 1:F
        % Reset n_iter
        n_iter = 1;
        
        
        % If concavity option on, converge until the right concavity is
        % reached
        while (range(f,1) <= range(f,2) - 30)
        %while (n_iter<2) || (convex_in~= 0 && convex_out(f,h)~=0 && convex_out(f,h) ~= convex_in && range(f,1) <= range(f,2) - 30) || ...
        %        (convex_in~= 0 && convex_out(f,h)~=0 && convex_out(f,h) == convex_in && test_res(n_iter) > test_res(n_iter-1) && range(f,1) <= range(f,2) - 30)
            %    (convex_in~= 0 && convex_out(f,h)~=0 && convex_out(f,h) == convex_in && range(f,1) <= range(f,2) - 30 && res(f,h) > .01)
            
            % Normalise
            xx_mat_norm(f,:) = normaliseFct(xx_mat(f,:), range(f,:));
            yy_mat_norm(f,:) = normaliseFct(yy_mat(f,:), range(f,:));
            
            % Fitting
            try
                range_fit = range(f,1):dx:range(f,2);
                [pfit(f,:,h),res(f,h)]=lsqcurvefit(@fun_exp_drift,p0(f,:),...
                    xx_mat_norm(f,range_fit),yy_mat_norm(f,range_fit),[],[],optsN);
            catch
                pfit(f,:,h) = p0(f,:);
                res(f,h) = 1e3;
                disp('Fitting error')
            end
            
            % Function recovery
            range_itv = range_dim(f,1,h):range_dim(f,2,h);
            out_exp(f,:,h) = denormaliseFct(real(fun_exp_drift(pfit(f,:,h),xx_mat_norm(f,:)')'),...
                range(f,:),yy_mat(f,:));            
            out_diff2 = diff(diff((out_exp(f,:,h))));
            convex_out(f,h) = sign(mean(out_diff2(:,range_itv(1):range_itv(end)-2)));
            test_convex(n_iter) = sign(mean(out_diff2(:,range_itv(1):range_itv(end)-2)));
            test_res(n_iter) = res(f,h);
            test_out_exp(n_iter,:) = out_exp(f,:,h);
            

            if convex_in == 0 % No condition on concavity, just finish
                break;
            else                
                %Plot 
%                 figure(15);
%                 subplot(2,1,1);
%                 plot(test_convex);
%                 subplot(2,1,2);
%                 plot(test_res);
%                 figclf(14); hold all;
%                 plot(time,out(f,:));
%                 plot(time,out_exp(f,:,h));
%                 plotValue(range(f,:),time,out_exp(f,:,h),'X')
            
                if (n_iter>1) && (test_convex(n_iter) == convex_in) && ...
                    (test_convex(n_iter-1) == convex_in) && ...
                    (test_res(n_iter)>test_res(n_iter-1))% End of the algo
                    %
                    out_exp(f,:,h) = test_out_exp(n_iter-1,:);
                    res(f,h) = test_res(n_iter-1);

                    break;
                end
                % Increment
                range(f,1) = range(f,1) + 30;
                n_iter = n_iter+1;
            end

%             % Plot - Troubleshoot
%             if n_iter>1 && (convex_in ~= 0 && convex_out(f,h) ~=0) && f>1000
%                 plot(time,out_exp(:,f,h));
%                 convex = convex_out(f,h)
%                 f
%             end
        end
        % Plot results
%         figclf(10); hold all;
%         plot(time,out(:,f),'-');
%         plot(time,out_exp(:,f,h),'--');
%         plot(xx_mat_norm(range_fit,f),yy_mat_norm(range_fit,f));
%         plot(xx_mat_norm(range_fit,f),fun_exp_drift(pfit(f,:,h),xx_mat_norm(range_fit,f)))
%         f    
    end
%     toc;
end

end

function [f]=fun_exp_drift(x,xdata)
    f = x(:,1).*(1-exp(-(xdata./x(:,2)).^x(:,3)));
%     f = x(:,1).*exp(x(:,2).*xdata) + x(:,3).*exp(x(:,4).*xdata);
end