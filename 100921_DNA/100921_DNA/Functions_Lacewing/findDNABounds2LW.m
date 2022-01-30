function [out_exp_i,out_exp_f,xi,xf] = findDNABounds2LW(out,time,range_in,x_stable,xm,x_tol,dx)
%findDNABounds find the bounds of the amplification for the sensor output
% Inputs:
% - out(NxP) : Sensor output array.
% - time(1xN) : Time vector;
% - range_in : Range for bounds. Range(1) is x_start and range(2) is larger
%   than xm.
% - xm(Px1) : Max derivative time stamp.
% - x_tol : Starting time stamp from data edge.
% - dx : Resolution of time stamp sweep.

if nargin<6
    x_tol = 20;
end
if nargin<7
    dx=10;
end

% Initialisation
N = size(out,2);
P = size(out,1);
Np = range_in(1,2); % Limited number of samples

% Parameters
tol = 1e-3;

% Start of the curve
% For each pixel
pixel_on_dna = zeros(P,1);
[out_exp_i,out_exp_f] = initialiseArray([P,N]);
[xi,xf] = initialiseArray([P,1]);
f = waitbar(0,'Drift compensation');
for p = 1:P
%     tic;
    % Interpolation range
    range_i = x_stable+x_tol:dx:xm(p);
    Xi = length(range_i);
    
%     % Debug
%     figclf(20); hold all;
%     plot(time,out(p,:),'linewidth',3);

    % Initialisation
    err = zeros(Xi,1);
    res = zeros(Xi,1);

    % Iteration of interpolation
    for x = 1:Xi
        range = [x_stable,range_i(x)];
        [out_exp_i(p,:),res(x)] = compensateExp2ParallelLW(out(p,:),time,range,[],10,tol);
        range_err = range_i(x)-5:range_i(x)+5;
        err(x) = immse(out(p,range_err),out_exp_i(p,range_err));
%         plot(time,out_exp_i(p,:));
    end

    % Find the dip in residue
%     figclf(21); 
%     subplot(2,1,1);plot(range_i,res)
%     subplot(2,1,2);plot(range_i,err)
    
    [dips_y,dips_x] = findpeaks(-err);
    if ~isempty(dips_x)
        pixel_on_dna(p) = 1;
        idx = dips_x(end);
    else
        idx = find(err<.3*err(end));
        idx = idx(end);
    end
    % Output final initial interpolation
    xi(p) = range_i(idx);
    range = [x_stable,range_i(idx)];
    [out_exp_i(p,:),~,pfit_i] = compensateExp2ParallelLW(out(p,:),time,range,[],10,1e-3);

%     % Debug - plot
%     figclf(22); hold all;
%     plot(time,out(p,:),'linewidth',3);
%     plot(time,out_exp_i(p,:))
    
    % Final interpolation
    if xm(p) < Np-10
        range_f = xm(p)+10:dx:Np;
        range = [range_f(1),range_f(end)];
        Xf = length(range_f);
        out_diff = diff((out_exp_i(p,x_stable:xm(p))));
        out_diff2 = diff(diff((out_exp_i(p,x_stable:xm(p)))));
        convex = sign(mean(out_diff2));
        [out_exp_f(p,:),~,~,~,range_out] = compensateExp2ParallelLW(out(p,:),time,range,pfit_i,10,tol,convex);
        xf(p) = range_out(1);

        % Find bounds
        % Linear aproximation around xm
        pol_m = polyfit(time(xm(p)-10:xm(p)+10),...
            out(p,xm(p)-10:xm(p)+10),1);
        y_pol_m = polyval(pol_m,time);
        %  Find deviation points
        [~,xi(p)] = min(abs(out_exp_i(p,:)-y_pol_m));
        [~,xf(p)] = min(abs(out_exp_f(p,:)-y_pol_m));
    else
        xi(p) = xm(p);
        xf(p) = xm(p);
    end
    
    waitbar(p/P,f)

end
close(f);

% figclf(22); hold all;
% plot(time,out(p,:),'linewidth',3);
% plot(time,out_exp_i(p,:),'linewidth',1);
% plot(time,out_exp_f(p,:),'linewidth',1);
% plotValue(xi,time,out(p,:),'X');
% plotValue(xf,time,out(p,:),'X');

% Find bounds
% Polyfit


end

