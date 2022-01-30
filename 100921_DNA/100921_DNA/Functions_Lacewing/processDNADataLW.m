function [time,out_smooth,ampl_curve_V,ampl_curve_pH,ampl_curve_lin,ampl_curve_sigm,...
          ampl_curve_norm,exp_comp,sigm_param,pixel_on,pixel_on_dna,xm,xi,xf,Ct,ampl,x_stable] ...
          = processDNADataLW(file_data,file_data_all,start_exp,end_exp,pH_init,options)
%extractDNAcurves extracts amplification curves for an array of sensors on
% the LACEWING platform.
% Input :
% - file_data : filename with average data and time stamps.
% - file_data_all : filename with all datasets.
% - start_exp : time range for loading time.
% - end_exp : time stamp for end of the reaction.
% - pH_init : initial pH of the reaction.
% - en_temp : 1 for detection of reaction start from the temperature
%   pixels, 0 for detection from the chemical pixels.
% - options : a cell of strings  to indicate options within the programme:
%       - 'startchem' : detection of reaction start from the chemical
%         pixels instead of temperature pixels.
%       - 'skipaverage' : do not run algo on average output.
%       - 'onlyaverage' : only run algo on average.
%       - 'wells' : number of wells wells.
%       - 'skipVrefcalib' : skip Vref calibration.
      
% Constants
[ROWS,COLS,P,coord_chem,ROWS_temp,COLS_temp,npt,coord_temp,coord_peltier] = getConstantsLW;
 
% Algorithm parameters
sens = 12;
minute = 60;
smooth_factor = 60;

% Option
if nargin<6
    options = {};
end
[pixel_spread,pixel_custom,W] = ...
    cellparam(options,{'spreadpixels','custompixels','wells'},[0,0,1]);

% Import
[time_all,out_all] = loadLW(file_data,file_data_all);

%% 1 - Calibration for active pixels
if ~cellContainsStr(options,'skipVrefcalib')
    pixel_on_Vref = compensateVrefLW(out_all);
else
    pixel_on_Vref  = 1:P;
end

% pixel_on_line = zeros(P,1);
% pixel_on_line(pixel_on_Vref) = 1;
% pixel_on_array = reshape(pixel_on_line,ROWS,COLS);
% figure(5); surf(pixel_on_array); view(2);

% Algorithm
if cellContainsStr(options,'raw')
    figclf(1);
    % Plot average of all pixels
    ax1 = subplot(4,4,[1 2]);
    plot(time_all/60,out_all);
    title('All pixels');
    setFigureOptions('','Voltage [mV]',12,'All pixels')
    ax2 = subplot(4,4,[5 6]);
    plot(time_all/60,out_all(coord_temp,:),'linewidth',2);
    linkaxes([ax1,ax2],'x');
    setFigureOptions('Time [min]','Voltage [mV]',12)
    return;
end

% Peltier voltage
[Vtemp,Vpeltier] = getPeltierLW(out_all,time_all,coord_peltier);
Vpeltier_reg = Vpeltier(199);

%% 2- Find the start of the experiments
if cellContainsStr(options,'startchem')
    en_temp = 0;
else
    en_temp = 1;
end
[range,x_stable,t_stable] = findStartLW(time_all,out_all,en_temp,start_exp,end_exp,W);
time = time_all(range(1):range(2))-time_all(range(1));
out = out_all(:,range(1):range(2));
N = size(out,2);
Vtemp = smooth(Vtemp(range(1):range(2)),smooth_factor);
Vpeltier = Vpeltier(range(1):range(2));

%% 3 - Signal pre-filtering
% 3a - Remove out of bound pixels
[pixel_on,pixel_off,Pon,Poff] = removeOutPixelsLW(out,pixel_on_Vref,x_stable);

% Coordinates of wells
thr = ROWS/W;
[coord_wells,wells] = divideArrayLW(thr);
for w = 1:W
    [~,coord_wells{w}] = intersect(pixel_on,coord_wells{w});
    Pon(w) = length(coord_wells{w});
end

% 3b/c - Output smoothing and subtracting with temperature output
out_smooth = zeros(P,N,W);
for w = 1:W
    for p = 1:Pon(w)
        out_smooth(pixel_on(coord_wells{w}(p)),:,w) = ...
            smooth(out(pixel_on(coord_wells{w}(p)),:),smooth_factor)-smooth(mean(out(coord_temp,:),1),smooth_factor) + 512;
        % Here add subtraction to closest temp pixel
    end
end

% 3d - Spatial averaging
% Initialisation
out_smooth_mov = out_smooth;
for w = 1:W
    % Moving average
    [~,out_smooth_mov(pixel_on(coord_wells{w}),:)] = movingAverageLW(out_smooth(pixel_on(coord_wells{w}),:),pixel_on(coord_wells{w}));
end
out_smooth = out_smooth_mov;

%% 4 - Spatial metrics

[xm,tm] = initialiseArray([P,1]);
for w = 1:W
    % 4a - Amplification peak 
    % Find max derivatives
    range = ones(Pon(w),1)*[x_stable N-60];
    [xm(pixel_on(coord_wells{w})),tm(pixel_on(coord_wells{w}))] = findDNAPeakLW(out_smooth(pixel_on(coord_wells{w}),:,w),time,4*smooth_factor,range);
    tm_mean(w) = mean(tm(pixel_on(coord_wells{w})));
%     % Cut dataset to twice average tm
%     Np = min(N,time_to_idx(time,tm_mean+1200));
%     % Remove tm outliers
    tm_total = tm(coord_wells{w})/60;
    tm_total = tm_total(tm_total>0);
%     pixel_on(idx) = [];
    Pon(w) = length(coord_wells{w});
end
% % Fit initial exponential
% range = ones(P,1)*time_to_idx(time,time(x_stable) + [60 mean(time(x0(pixel_on)))/2]);
% out_exp_init = zeros(P,N);
% [out_exp_init(pixel_on,:)] = compensateExp2ParallelLW(out_smooth(pixel_on,:),time,range);
% out_comp = out_smooth - out_exp_init;
% % Find accurate maximum of derivative
% [xm,tm] = initialiseArray([P,1]);
% [xm(pixel_on),tm(pixel_on)] = findDNAPeakLW(out_comp(pixel_on,:),time,smooth_factor,x_stable,[time(x_stable) time(end)-3*60]);

figclf(1);
subplot(4,4,[1 5]);
surfArrayLW(tm,'reshape');
setFigureOptions('Pixel X','Pixel Y',12,'Max derivative')
subplot(4,4,3);
hist(tm_total,100);
setFigureOptions('','',12,'Max derivative')
xlabel('Time [min]');
ylabel('Occurences');

% % Plot all max derivatives
% figclf(1); subplot(4,4,[1 5]);
% surfArrayLW(tm,'reshape');
% setFigureOptions('Pixel X','Pixel Y',12,'Max derivative')
% colormap parula;

% 4b - Exponential correlation
res_exp_log = zeros(P,1);
for w = 1:W
    % Full interpolation
    range = ones(P,1)*[x_stable N];
    res_exp = zeros(P,1);
    [~,res_exp(pixel_on(coord_wells{w}))] = compensateExp2ParallelLW(out_smooth(pixel_on(coord_wells{w}),:,w),time,range,[],30,1e-3);
    res_exp_log(pixel_on(coord_wells{w})) = -log10(res_exp(pixel_on(coord_wells{w})));
end

figure(1);
subplot(4,4,[2 6]);
surfArrayLW(res_exp_log,'reshape'); 
setFigureOptions('Pixel X','Pixel Y',12,'Residual')
subplot(4,4,7);
hist(res_exp_log(pixel_on),100);
setFigureOptions('','',12,'Residual')


%% 5 - Clustering for each well

% % Cluster
% coord_cluster = cell(W,1);
% for w = 1:W
%     [coord_cluster{w},~,C] = ...
%         clusterLW(log10(res_full(coord_wells{w})),pixel_on(coord_wells{w}),10);
% %     [~,idx] = max(cellfun('size',coord_cluster{w},1));
%     for c = 1:C
%         coord_cluster{w}{c} = coord_wells{w}(coord_cluster{w}{c});
%     end
% end

for w = 1:W
    % Cluster based on residual
    [coord_cluster{w},~,C] = clusterLW(res_exp_log(pixel_on(coord_wells{w})),pixel_on(coord_wells{w}),10);
    % Note : build up for several wells
    % Decision on 'active' and 'barrier' clusters
    [coord_cluster_final{w},res_cluster_final{w},tm_cluster_final{w}, ...
        density_cluster_final{w}, res_cluster{w}, tm_cluster{w}, density_cluster{w}] =...
    classifyClusterLW(coord_cluster{w},time,pixel_on(coord_wells{w}),res_exp_log,tm,x_stable);
    % DBSCAN on remaining active cluster
    pixel_on_cluster{w} = findWellLW(coord_cluster_final{w},1,pixel_on(coord_wells{w}));
end

figure(1);
subplot(4,4,[4 8]);
array=zeros(P,1);
for w = 1:W
    array(pixel_on_cluster{w})=w;
end
surfArrayLW(array,'reshape')
setFigureOptions('Pixel X','Pixel Y',12,'On pixels')

% tm_mean = tm_cluster_final(1)*60;
% Np = min(N,time_to_idx(time,tm_mean+1200));

% figclf(22);
% array = zeros(P,1);
% for c = 1:C
%     array(pixel_on(coord_cluster{c}))=c;
% end
% surfArrayLW(array,'reshape')
% 
% Plot average curves for each cluster
% figclf(26);
% for c = 1:C
%     subplot(C,1,c);
%     plot(time/60,mean(out_smooth(pixel_on(coord_cluster{c}),:)),'Linewidth',2);
%     setFigureOptions('Time [min]','Output [mV]',16,...
%         ['P = ' num2str(length(coord_cluster{c})) ...
%         ',residual = ' num2str(res_cluster(c)) ...
%         ', tm = ' num2str(tm_cluster(c),3)...
%         ', density = ' num2str(density_cluster(c))])
% end

% 

% 
% % Sort positives and negatives - 0 unstable, 1 positive, 2 negative
% pixel_dna = zeros(size(pixel_on));
% for w = 1:W
%     pixel_dna(coord_cluster{w}) = 1;
%     idx = union(find(tm(coord_cluster{w})<range(1)+2*60),find(tm(coord_cluster{w})>range(2)-2*60));
%     pixel_dna(coord_cluster{w}(idx)) = 2;
% end
% 
% % Plot DNA heat map
% figure(1);
% subplot(4,4,[2 6]);
% surfArrayLW(pixel_dna,'fill',pixel_on);
% setFigureOptions('Pixel X','Pixel Y',12,'Post-clustering')

% Average of array
for w = 1:W
    out_smooth_av{w} = mean(out_smooth(pixel_on_cluster{w},:,w),1);%mean(out_smooth(:,coord_wells{1}),2);
end

% if w>1
%     out_sub = diff(out_smooth_av,1,2); % Subtraction
%     out_sub = out_sub-out_sub(1)+500; % Levelling 
% else
%     out_sub = 0;
% end

% Plot curves
figure(1);
for w = 1:W
    ax(1) = subplot(4,4,9+w-1); hold all;
    plot(time/60,out_smooth_av{w},'Linewidth',2);
    setFigureOptions('','Voltage [mV]',12,'Chemical voltage',[0 time(end)/60]);
%     xlim([0 time(end)/60]);
%     legend('Control','Sample');
    % ax(2) = subplot(6,4,17);
    % plot(time/60,out_sub,'Linewidth',2);
    % setFigureOptions('','Voltage [mV]',12,'Subtraction')
    ax(3) = subplot(4,4,13+w-1);
    plot(time/60,smooth(Vtemp,smooth_factor),'linewidth',2);
    linkaxes(ax,'x');
    setFigureOptions('Time [min]','Voltage [mV]',12,'Temperature voltage')
end
% Cut experiment to double xm
% out_smooth_av = out_smooth_av(1:min(length(time),2*x0));
% time = time(1:min(length(time),2*x0));

if ~cellContainsStr(options,'skipaverage')
    for w = 1:W
        % Extract curves
        range = [x_stable,N];
        [xm_av{w},xi_av{w},xf_av{w},ampl_curve_V_av{w},exp_comp_av{w},ampl_curve_pH_av{w},...
            ampl_curve_lin_av{w},ampl_curve_sigm_av{w},ampl_curve_norm{w},time_sigm{w},outcome{w},sigm_param_av{w},range] ...
            = extractDNACurvesLW(out_smooth_av{w},time,smooth_factor,range,pH_init,0);

        figure(1); 
        subplot(4,4,9+w-1);
        plot(time/60,exp_comp_av{w},'--','Linewidth',2);
        setFigureOptions('Time [min]','Voltage [mV]',12,'Chemical voltage',[0 time(end)/60]);
        legend('Sensor output','Interpolated sensor drift','Location','SouthEast')

        subplot(4,4,10); hold all;
        plot(time/60,ampl_curve_V_av{w},'Linewidth',2);
        setFigureOptions('Time [min]','Voltage [mV]',12,'Compensated voltage output',[0 time(end)/60])

        subplot(4,4,11); hold all;
        plot(time/60,ampl_curve_lin_av{w},'linewidth',2);
        h=plot(time_sigm{w}/60,ampl_curve_sigm_av{w},'--','linewidth',1.5);
        setFigureOptions('Time [min]','Proton count',12,'Amplification',[0 time_sigm{w}(end)/60])
        legend('Amplification curve','Sigmoid fitting','Location','NorthWest')

        subplot(4,4,15); hold all;
        plot(time_sigm{w}/60,ampl_curve_norm{w},'Linewidth',2);
        setFigureOptions('Time [min]','Voltage [mV]',12,'Normalised amplification',[0 time_sigm{w}(end)/60])

        [~,x_Ct] = min(abs(ampl_curve_norm{w}-.1));
        Ct(w) = time_sigm{w}(x_Ct)/60; % CHANGE TO TIME_SIGM
    end
    
%     for w = 1:W
%         if outcome(w) == 1
%             idx = find(ampl_curve_sigm_av(:,w)>ampl_curve_sigm_av(1,w)+.2*(ampl_curve_sigm_av(end,w)-ampl_curve_sigm_av(1,w)));
%             disp(['TTP for well ' num2str(w) ' is ' num2str(time(idx(1))/60) ' min'])
%             legendStr{w}= ['TTP for well ' num2str(w) ' = ' num2str(time(idx(1))/60) ' min'];
%         elseif outcome(w) == 2
%             disp(['Well ' num2str(w) ' is negative']);
%             legendStr{w}=['Well ' num2str(w) ' is negative'];
%         elseif outcome(w) == 0
%             disp(['Well ' num2str(w) ' is unstable']);
%             legendStr{w}=['Well ' num2str(w) ' is unstable'];
%         end
%     end
%     legend(legendStr);
%     
%     % Subtraction
%     if w>1
%         [xm_sub,xi_sub,xf_sub,ampl_curve_V_sub,exp_comp_sub,~,...
%             ampl_curve_lin_sub,ampl_curve_sigm_sub,outcome,sigm_param_sub] ...
%             = extractDNACurvesLW(out_sub,time,smooth_factor,x_stable,pH_init);
%     else
%         [xm_sub,xi_sub,xf_sub] = initialiseArray([1,1]);
%         [ampl_curve_V_sub,exp_comp_sub,...
%             ampl_curve_lin_sub,ampl_curve_sigm_sub,outcome,sigm_param_sub] = initialiseArray(size(time));
%     end
%     figure(1); subplot(4,4,15); hold all;
%     plot(time/60,ampl_curve_lin_sub,'linewidth',2);
%     plot(time/60,ampl_curve_sigm_sub,'--','linewidth',1.5);
%     setFigureOptions('Time [min]','Normalised proton count',12,'Amplification after subtraction')
%     if outcome == 1
%         idx = find(ampl_curve_lin_sub(:)>ampl_curve_lin_sub(1)+.2*(ampl_curve_lin_sub(end)-ampl_curve_lin_sub(1)));
%         disp(['TTP for subtracted curve is ' num2str(time(idx(1))/60) ' min'])
%     elseif outcome == 2
%         disp(['Subtraction is negative']);
%     elseif outcome == 0
%         disp(['Subtraction is unstable']);
%     end

    % Results
    for w = 1:W
        disp(['Well ' num2str(w) ' residual = ' num2str(mean(res_exp_log(pixel_on_cluster{w})),3)]);
        disp(['Well ' num2str(w) ' Ct = ' num2str(Ct(w),3) 'min']);
    end

end



if cellContainsStr(options,'onlyaverage')
    out_smooth=out_smooth_av;
    ampl_curve_V=ampl_curve_V_av;
    ampl_curve_pH=ampl_curve_pH_av;
    ampl_curve_lin=ampl_curve_lin_av;
    ampl_curve_sigm=ampl_curve_sigm_av;
    exp_comp=exp_comp_av;
    sigm_param=sigm_param_av;
    pixel_on_dna=0;
    xm=xm_av;
    xi=xi_av;
    xf=xf_av;
    ampl=0;
    return;
end

% Spatial
% Initialisation
[xm,xi,xf] = initialiseArray([P,1]);
[ampl_curve_V,exp_comp,ampl_curve_pH,ampl_curve_lin,ampl_curve_sigm,...
    ampl_curve_norm] = initialiseArray([P,N]);

% Extract curves
range = [x_stable,N];
[xm(pixel_on_cluster{1}),xi(pixel_on_cluster{1}),xf(pixel_on_cluster{1}),...
    ampl_curve_V(pixel_on_cluster{1},:),exp_comp(pixel_on_cluster{1},:),...
    ampl_curve_pH(pixel_on_cluster{1},:),ampl_curve_lin(pixel_on_cluster{1},:),...
    ampl_curve_sigm(pixel_on_cluster{1},:),ampl_curve_norm(pixel_on_cluster{1},:),...
    time_sigm,~,sigm_param,range(pixel_on_cluster{1},:)] ...
    = extractDNACurvesLW(out_smooth(pixel_on_cluster{1},:),time,smooth_factor,range,pH_init,0);

figclf(10);
plot(time,ampl_curve_norm)
size('Dengue 10^6 - All amplification curves')
title('Dengue 10^6 - All amplification curves')
set(gca,'fontsize',14)
xlabel('Time [s]')
ylabel('Intensity [A.U.]')

figure(11); surf(reshape(tm,ROWS,COLS))
view(2);
tm(pixel_on_cluster{1})
mean(tm(pixel_on_cluster{1}))
std(tm(pixel_on_cluster{1}))
caxis([mean(tm(pixel_on_cluster{1}))-2*std(tm(pixel_on_cluster{1}))
       mean(tm(pixel_on_cluster{1}))+2*std(tm(pixel_on_cluster{1}))]);
colorbar

% % Post-filtering
% pixel_dna(sigm_param(1,:)<4) = 2;
% pixel_dna(sigm_param(3,:)>1) = 2;
% pixel_dna(ampl_curve_sigm(end,:)-ampl_curve_sigm(1,:)>1e-7) = 2;
% pixel_off = union(find(sigm_param(1,:)<4),find(sigm_param(3,:)>1));
% pixel_off = union(pixel_off,find(ampl_curve_sigm(end,:)-ampl_curve_sigm(1,:)>1e-7));
% pixel_on_dna(pixel_off) = [];
% sigm_param(:,pixel_off) = [];

% Normalise
ampl_curve_norm = (ampl_curve_sigm-ampl_curve_sigm(1,:))./(ampl_curve_sigm(end,:)-ampl_curve_sigm(1,:));
% ampl_curve_norm(:,pixel_off) = [];

figure(1);
subplot(4,4,[11 12 15 16]); hold off;
plot(time/60,ampl_curve_norm,'linewidth',2);
setFigureOptions('Time [min]','Proton concentration [M]',12,'Sigmoidal fitting');

% Default value
pixel_on_dna=pixel_on;
ampl=0;
Ct = 0;

% % Plot DNA heat map
% figure(1);
% subplot(4,4,[2 6]);
% surfArrayLW(pixel_dna,'fill',pixel_on);
% setFigureOptions('Pixel X','Pixel Y',12,'Max derivative')
% colormap jet;

% % Plot Ct
% figure(1);
% subplot(4,4,[3 7]);
% surfArrayLW(sigm_param(1,:),'fill',pixel_on(pixel_on_dna));
% setFigureOptions('Pixel X','Pixel Y',12,'Max derivative')
% colormap jet;

% 
% 
% % Feature extraction
% colormap jet;
% feature_label = {'Ct','Slope','Amplitude','Asymmetry'};
% figclf(2);
% for i = 1:4
%     subplot(2,2,i);
%     surfArrayLW(sigm_param(i,pixel_dna==1),'fill',pixel_on(pixel_dna==1));
%     setFigureOptions('Pixel X','Pixel Y',12,feature_label{i})
%     view(2);
%     colorbar;
% end

% Average amplitude and TTP
% for w = 1:W
%     [~,coord_wells_dna{w}] = intersect(pixel_on_dna,coord_wells{w});
%     ampl(w) = mean(sigm_param(3,coord_wells_dna{w}))
%     Ct(w) = mean(sigm_param(1,coord_wells_dna{w}));
%     disp(['Average pH change in well ' num2str(w) ' is ' num2str(ampl(w))])
%     disp(['Average Ct in well ' num2str(w) ' is ' num2str(Ct(w)) ' mins'])
% end

disp('Over')
end