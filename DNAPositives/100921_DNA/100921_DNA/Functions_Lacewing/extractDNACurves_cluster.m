function [] = extractDNACurves_cluster(file_data,file_data_all,start_exp,end_exp,pH_init,t0,en_temp,pixel_on)

% Constants
[ROWS,COLS,np,coord_chem,ROWS_temp,COLS_temp,npt,coord_temp,coord_peltier] = getConstantsLW;

% Algorithm option
smooth_factor = 40;
sens = 12;
minute = 60;

% Import
[time_all,out_all] = loadLW(file_data,file_data_all);

%% Calibration for active pixels
pixel_on_Vref = compensateVrefLW(out_all);
pixel_on_Vref = 1:np;

%% Algorithm
figclf(1);
% Plot average of all pixels
ax1 = subplot(4,4,[1 2]);
plot(time_all/60,out_all);
title('All pixels');
setFigureOptions('','Voltage [mV]',12,'All pixels')
ax2 = subplot(4,4,[5 6]);
plot(time_all/60,out_all(:,coord_temp),'linewidth',2);
linkaxes([ax1,ax2],'x');
setFigureOptions('Time [min]','Voltage [mV]',12)

% Peltier voltage
[Vtemp,Vpeltier] = getPeltierLW(out_all,time_all,coord_peltier);
Vpeltier_reg = Vpeltier(199);

% Find the start of the experiments
[range,x_stable,t_stable] = findStartLW(time_all,out_all,en_temp,start_exp,end_exp);
time = time_all(range(1):range(2))-time_all(range(1));
out = out_all(range(1):range(2),:);
n = size(out,1);
Vtemp = smooth(Vtemp(range(1):range(2)),smooth_factor);
Vpeltier = Vpeltier(range(1):range(2));

% figure(2); clf(2);
% ax1=subplot(2,1,1);
% plot(time,Vtemp);
% ax2=subplot(2,1,2);
% plot(time,Vpeltier);
% linkaxes([ax1,ax2],'x');
% setFigureOptions('Time [min]','Voltage [mV]',12,'Temperature pixels')

% Remove out of bound pixels
% [pixel_on, pixel_off, npon,npoff] = removeOutPixelsLW(out,pixel_on_Vref,x_stable);
pixel_off = setdiff(1:np,pixel_on);
npon = length(pixel_on);
npoff = length(pixel_off);

% Output smoothing and subtracting with temperature output
out_smooth = zeros(n,npon);
for p = 1:npon
    out_smooth(:,p) = ...
        smooth(out(:,pixel_on(p)),smooth_factor)-smooth(mean(out(:,coord_temp),2),smooth_factor) + 512;
    % Here add subtraction to closest temp pixel
end

% Spatial averaging
out_smooth_mov = zeros(size(out_smooth));
for i = 1:n
    [~,out_smooth_mov(i,:)] = movingAverageLW(out_smooth(i,:),pixel_on);
end
out_smooth = out_smooth_mov;

% Pixel on for DNA detection
pixel_on_dna = pixel_on;

% Fit the first three minutes to an exponential and subtract
out_smooth_av = mean(out_smooth,2);
% [out_exp,out_comp_av,~] = compensateExp2LW(out_smooth_av,time,[time(x_stable),3*60],1,2);

% Plot average of all pixels
figclf(1);
ax1 = subplot(4,4,[1 2]);
plot(time/60,out_smooth_av);
title('All pixels');
setFigureOptions('','Voltage [mV]',12,'All pixels')
ax2 = subplot(4,4,[5 6]);
plot(time/60,smooth(Vtemp,smooth_factor),'linewidth',2);
linkaxes([ax1,ax2],'x');
setFigureOptions('Time [min]','Voltage [mV]',12)

% Find maximum of derivative
x0 = findDNAPeakLW(out_smooth_av,time,smooth_factor,x_stable,[time(x_stable) time(end)-3*60]);
disp(['Mean derivative maximum at ' num2str(time(x0)/60) ' min']);

% Cut experiment to double xm
% out_smooth_av = out_smooth_av(1:min(length(time),2*x0));
% time = time(1:min(length(time),2*x0));

% Plot
figure(1); 
subplot(4,4,[3 7]); hold all;
plot(time/60,out_smooth_av,'linewidth',2);
plot(time(x0)/60*ones(1,100),linspace(out_smooth_av(x0)-10,out_smooth_av(x0)+10,100),'k');
setFigureOptions('Time [min]','Voltage [mV]',12,'Average output')

% Comment average
% % Find amplification region for average
% [xi_av,xf_av,out_exp_i,out_exp_f] = findDNABoundsLW(out_smooth_av,time,x_stable,x0,smooth_factor);
% % Plot
% figure(1); subplot(4,4,[3 7]); hold all;
% plot(time(xi_av)/60*ones(1,100),linspace(out_smooth_av(xi_av)-10,out_smooth_av(xi_av)+10,100));
% plot(time(xf_av)/60*ones(1,100),linspace(out_smooth_av(xf_av)-10,out_smooth_av(xf_av)+10,100));
% plot(time(x0)/60*ones(1,100),linspace(out_smooth_av(x0)-10,out_smooth_av(x0)+10,100));
% % plot(time/60,out_exp_i,'--','linewidth',2);
% % plot(time/60,out_exp_f,'--','linewidth',2);
% 
% % Drift compensation
% [ampl_curve_V_av,exp_comp_av] = compensateDriftLW(time,out_smooth_av,xi_av,xf_av,out_exp_i,out_exp_f,x_stable);
% figure(1); 
% subplot(4,4,[3 7]);
% plot(time/60,exp_comp_av,'--','linewidth',2);
% subplot(4,4,4); hold all;
% plot(time/60,ampl_curve_V_av,'linewidth',2);
% plot(time(x0)/60*ones(1,100),linspace(ampl_curve_V_av(x0)-2,ampl_curve_V_av(x0)+2,100),'k');
% plot(time(xi_av)/60*ones(1,100),linspace(ampl_curve_V_av(xi_av)-2,ampl_curve_V_av(xi_av)+2,100));
% plot(time(xf_av)/60*ones(1,100),linspace(ampl_curve_V_av(xf_av)-2,ampl_curve_V_av(xf_av)+2,100));
% setFigureOptions('Time [min]','Voltage [mV]',12,'Compensated output')
% 
% [ampl_curve_lin_av] = linAmplCurveLW(ampl_curve_V_av,x_stable,pH_init,smooth_factor);
% 
% % figclf(3);
% % for i = 1:length(pH_init)
% %     subplot(4,4,i);
% %     plot(time,ampl_curve_lin_av(i,:),'linewidth',2);
% %     title(['pH =' num2str(pH_init(i))])
% % end
%     
% subplot(4,4,8); hold all;
% plot(time/60,ampl_curve_lin_av,'linewidth',2);
% setFigureOptions('Time [min]','Voltage [mV]',12,'Compensated output')
%End comment

%% Then for each pixel

% Initialisation
[out_comp] = initialiseArray([n,np]);
[rsquare,xm,tm] = initialiseArray([np,1]);
% Constants
npon = length(pixel_on_dna);
[xm tm] = findDNAPeakLW(out_smooth,time,4*smooth_factor,x_stable,[time(1),time(end)-300]);

% Uncomment to plot;
figure(1);
subplot(4,4,[9 13]);
surfArrayLW(tm/60,'fill',pixel_on_dna);
setFigureOptions('Pixel X','Pixel Y',12,'Max derivative')

% % Consider no amplification if xm is out of bound
% idx = union(find(tm(pixel_on_dna) < 60*3),find(tm(pixel_on_dna) > time(end)-3*60));
% if ~isempty(idx)
%     xm(pixel_on_dna(idx)) = 0;
%     tm(pixel_on_dna(idx)) = 0;
%     pixel_on_dna(idx) = [];
% end

% Initialisation
[out_exp_i,out_exp_f,ampl_curve_V,ampl_curve_pH,ampl_curve_lin,exp_comp,ampl_curve_sigm] = initialiseArray([n,npon]);
[xi,xf,Ct] = initialiseArray([npon,1]);
% Waitbar
f = waitbar(0,'Processing amplification curves');
for p = 1:npon
    % Only process if the amplification peak looks alright
    if xm(p) > 5
        % Find amplification region for average
        [xi(p),xf(p),out_exp_i(:,p),out_exp_f(:,p)] = ...
            findDNABoundsLW(out_smooth(:,p),time,x_stable,x0,smooth_factor);
        
        % Drift compensation
        [ampl_curve_V(:,p),exp_comp(:,p)] = ...
            compensateDriftLW(time,out_smooth(:,p),xi(p),xf(p),...
            out_exp_i(:,p),out_exp_f(:,p),x_stable);

        % Linearisation
        [ampl_curve_lin(:,p),ampl_curve_pH(:,p)] = linAmplCurveLW(ampl_curve_V(:,p),x_stable,pH_init,smooth_factor);

        % Sigmoidal fitting
        [ampl_curve_sigm(:,p),sigmoid_fit,sigmoid_gof] = fitSigmoidLW(ampl_curve_lin(:,p),time,x_stable);

        thr = .1;
        thr_curve = ampl_curve_sigm(1,p)+.1*(ampl_curve_sigm(end,p)-ampl_curve_sigm(1,p));
        idx = find(ampl_curve_sigm(:,p)>thr_curve);
        if ~isempty(idx)
            Ct(p) = time(idx(1));
        else
            Ct(p) = 0;
        end
    end
    
    % Sigmoidal feature extraction
%     sigm_param = extractSigmParamLW(sigmoid_fit,sigmoid_gof,time);
    
%     figclf(6);
%     subplot(3,1,1); hold all;
%     plot(time,out_smooth(:,p));
%     plot(time,exp_comp(:,p));
%     plotValue(xi(p),time,out_smooth(:,p),'X','--');
%     plotValue(xm(p),time,out_smooth(:,p),'X','--');
%     plotValue(xf(p),time,out_smooth(:,p),'X','--');
%     
%     subplot(3,1,2); hold all;
%     plot(time,ampl_curve_V(:,p));
%     plotValue(xi(p),time,ampl_curve_V(:,p),'X','--');
%     plotValue(xm(p),time,ampl_curve_V(:,p),'X','--');
%     plotValue(xf(p),time,ampl_curve_V(:,p),'X','--');
% 
%     subplot(3,1,3); hold all;
%     plot(time,ampl_curve_lin(:,p));
%     plotValue(xi(p),time,ampl_curve_lin(:,p),'X','--');
%     plotValue(xm(p),time,ampl_curve_lin(:,p),'X','--');
%     plotValue(xf(p),time,ampl_curve_lin(:,p),'X','--');
    
%     waitforbuttonpress;
    waitbar(p/npon, f,sprintf('Pixel number %f / %f',p,npon));
end
close(f);

figure(1);
subplot(4,4,[10 11 14 15]);
plot(time/60,ampl_curve_sigm);
setFigureOptions('Pixel X','Pixel Y',12,'Amplification curve')

subplot(4,4,[12 16]);
surf(idx_to_array(Ct,pixel_on,ROWS,COLS));
view(2);
colorbar;

disp('Over')
end