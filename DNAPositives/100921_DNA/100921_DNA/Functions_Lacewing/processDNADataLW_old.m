function [time,out_smooth,ampl_curve_V,ampl_curve_pH,ampl_curve_lin,ampl_curve_sigm,...
          exp_comp,sigm_param,pixel_on,pixel_on_dna,xm,xi,xf,Ct,ampl] ...
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
smooth_factor = 40;

% Option
if nargin<6
    options = {};
end
[pixel_spread,pixel_custom,nw] = ...
    cellparam(options,{'spreadpixels','custompixels','wells'},[0,0,1]);

% Import
[time_all,out_all] = loadLW(file_data,file_data_all);

%% 1 - Calibration for active pixels
if ~cellContainsStr(options,'skipVrefcalib')
    pixel_on_Vref = compensateVrefLW(out_all);
else
    pixel_on_Vref  = 1:P;
end

% Algorithm
if cellContainsStr(options,'raw')
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
[range,x_stable,t_stable] = findStartLW(time_all,out_all,en_temp,start_exp,end_exp,nw);
time = time_all(range(1):range(2))-time_all(range(1));
out = out_all(range(1):range(2),:);
n = size(out,1);
Vtemp = smooth(Vtemp(range(1):range(2)),smooth_factor);
Vpeltier = Vpeltier(range(1):range(2));

%% 3 - Signal pre-filtering
% 3a - Remove out of bound pixels
[pixel_on, pixel_off, npon,npoff] = removeOutPixelsLW(out,pixel_on_Vref,x_stable);

% Coordinates of wells
thr = ROWS/nw;
[coord_wells,wells] = divideArrayLW(thr);
for w = 1:nw
    [~,coord_wells{w}] = intersect(pixel_on,coord_wells{w});
end

% 3b/c - Output smoothing and subtracting with temperature output
out_smooth = zeros(n,npon);
for p = 1:npon
    out_smooth(:,p) = ...
        smooth(out(:,pixel_on(p)),smooth_factor)-smooth(mean(out(:,coord_temp),2),smooth_factor) + 512;
    % Here add subtraction to closest temp pixel
end

% 3d - Spatial averaging
out_smooth_mov = zeros(size(out_smooth));
for i = 1:n
    [~,out_smooth_mov(i,:)] = movingAverageLW(out_smooth(i,:),pixel_on);
end
out_smooth = out_smooth_mov;



%% 4 - Clustering for each well
nmin = 5; % Min number of pixels in a cluster
coord_cluster = cell(nw,1);
for w = 1:nw
    coord_cluster{w} = ...
        findRegionLW(tm(coord_wells{w}),pixel_on(coord_wells{w}),thr,nmin,10,{'kmeans'});
    [~,idx] = max(cellfun('size',coord_cluster{w},1));
    coord_cluster{w} = coord_wells{w}(coord_cluster{w}{idx});
end

% Sort positives and negatives - 0 unstable, 1 positive, 2 negative
pixel_dna = zeros(size(pixel_on));
for w = 1:nw
    pixel_dna(coord_cluster{w}) = 1;
    idx = union(find(tm(coord_cluster{w})<range(1)+2*60),find(tm(coord_cluster{w})>range(2)-2*60));
    pixel_dna(coord_cluster{w}(idx)) = 2;
end

% Plot DNA heat map
figure(1);
subplot(4,4,[2 6]);
surfArrayLW(pixel_dna,'fill',pixel_on);
setFigureOptions('Pixel X','Pixel Y',12,'Post-clustering')

% Average of array
out_smooth_av = zeros(length(time),nw);
for w = 1:nw
    out_smooth_av(:,w) = mean(out_smooth(:,coord_cluster{w}),2);%mean(out_smooth(:,coord_wells{1}),2);
end
if w>1
    out_sub = diff(out_smooth_av,1,2); % Subtraction
    out_sub = out_sub-out_sub(1)+500; % Levelling 
else
    out_sub = 0;
end

% Plot curves
figure(1);
ax(1) = subplot(6,4,13);
plot(time/60,out_smooth_av,'Linewidth',2);
setFigureOptions('','Voltage [mV]',12,'Clusters');
legend('Control','Sample');
ax(2) = subplot(6,4,17);
plot(time/60,out_sub,'Linewidth',2);
setFigureOptions('','Voltage [mV]',12,'Subtraction')
ax(3) = subplot(6,4,21);
plot(time/60,smooth(Vtemp,smooth_factor),'linewidth',2);
linkaxes(ax,'x');
setFigureOptions('Time [min]','Voltage [mV]',12,'Temperature')

% Cut experiment to double xm
% out_smooth_av = out_smooth_av(1:min(length(time),2*x0));
% time = time(1:min(length(time),2*x0));

%% 4 - Find derivative

% Constants
npon = length(pixel_on);

% Find max derivatives
range = [time(x_stable) time(end)-3*60];
[x0,t0] = findDNAPeakLW(out_smooth,time,4*smooth_factor,x_stable,range);
% Fit initial exponential
range = ones(P,1)*time_to_idx(time,time(x_stable) + [60 mean(time(x0))/2]);
[out_exp_init,res_f] = compensateExp2ParallelLW(out_smooth,time,range);
out_comp = out_smooth - out_exp_init;
% Find accurate maximum of derivative
[xm,tm] = findDNAPeakLW(out_comp,time,smooth_factor,x_stable,[time(x_stable) time(end)-3*60]);

% Plot all max derivatives
figclf(1); subplot(4,4,[1 5]);
surfArrayLW(tm/60,'fill',pixel_on);
setFigureOptions('Pixel X','Pixel Y',12,'Max derivative')
colormap parula;

if ~cellContainsStr(options,'skipaverage')
    % Extract curves
    [xm_av,xi_av,xf_av,ampl_curve_V_av,exp_comp_av,ampl_curve_pH_av,...
        ampl_curve_lin_av,ampl_curve_sigm_av,outcome,sigm_param_av] ...
        = extractDNACurvesLW(out_smooth_av,time,smooth_factor,x_stable,pH_init,0);
    
    figure(1); subplot(4,4,10); hold all;
    plot(time,out_smooth_av,'Linewidth',2);
    plot(time,exp_comp_av,'--');
    setFigureOptions('','Voltage [mV]',12,'')
    
    figure(1); subplot(4,4,14); hold all;
    plot(time,ampl_curve_V_av,'Linewidth',2);
    setFigureOptions('Time [min]','Voltage [mV]',12,'')
    
    figure(1); subplot(4,4,11); hold all;
    plot(time/60,ampl_curve_lin_av,'linewidth',2);
    h=plot(time/60,ampl_curve_sigm_av,'--','linewidth',1.5);
    setFigureOptions('','Normalised proton count',12,'Amplification in each well')
    for w = 1:nw
        if outcome(w) == 1
            idx = find(ampl_curve_sigm_av(:,w)>ampl_curve_sigm_av(1,w)+.2*(ampl_curve_sigm_av(end,w)-ampl_curve_sigm_av(1,w)));
            disp(['TTP for well ' num2str(w) ' is ' num2str(time(idx(1))/60) ' min'])
            legendStr{w}= ['TTP for well ' num2str(w) ' = ' num2str(time(idx(1))/60) ' min'];
        elseif outcome(w) == 2
            disp(['Well ' num2str(w) ' is negative']);
            legendStr{w}=['Well ' num2str(w) ' is negative'];
        elseif outcome(w) == 0
            disp(['Well ' num2str(w) ' is unstable']);
            legendStr{w}=['Well ' num2str(w) ' is unstable'];
        end
    end
    legend(legendStr);
    
    % Subtraction
    if w>1
        [xm_sub,xi_sub,xf_sub,ampl_curve_V_sub,exp_comp_sub,~,...
            ampl_curve_lin_sub,ampl_curve_sigm_sub,outcome,sigm_param_sub] ...
            = extractDNACurvesLW(out_sub,time,smooth_factor,x_stable,pH_init);
    else
        [xm_sub,xi_sub,xf_sub] = initialiseArray([1,1]);
        [ampl_curve_V_sub,exp_comp_sub,...
            ampl_curve_lin_sub,ampl_curve_sigm_sub,outcome,sigm_param_sub] = initialiseArray(size(time));
    end
    figure(1); subplot(4,4,15); hold all;
    plot(time/60,ampl_curve_lin_sub,'linewidth',2);
    plot(time/60,ampl_curve_sigm_sub,'--','linewidth',1.5);
    setFigureOptions('Time [min]','Normalised proton count',12,'Amplification after subtraction')
    if outcome == 1
        idx = find(ampl_curve_lin_sub(:)>ampl_curve_lin_sub(1)+.2*(ampl_curve_lin_sub(end)-ampl_curve_lin_sub(1)));
        disp(['TTP for subtracted curve is ' num2str(time(idx(1))/60) ' min'])
    elseif outcome == 2
        disp(['Subtraction is negative']);
    elseif outcome == 0
        disp(['Subtraction is unstable']);
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
    Ct = num2str(time(idx(1))/60);
    ampl=0;
    return;
end

% Spatial
pixel_on_dna = find(pixel_dna==1);
[xm,xi,xf,ampl_curve_V,exp_comp,ampl_curve_pH,ampl_curve_lin,ampl_curve_sigm,outcome,sigm_param] ...
    = extractDNACurvesLW(out_smooth(:,pixel_on_dna),time,smooth_factor,x_stable,pH_init,1);

% % Post-filtering
pixel_dna(sigm_param(1,:)<4) = 2;
pixel_dna(sigm_param(3,:)>1) = 2;
pixel_dna(ampl_curve_sigm(end,:)-ampl_curve_sigm(1,:)>1e-7) = 2;
pixel_off = union(find(sigm_param(1,:)<4),find(sigm_param(3,:)>1));
pixel_off = union(pixel_off,find(ampl_curve_sigm(end,:)-ampl_curve_sigm(1,:)>1e-7));
pixel_on_dna(pixel_off) = [];
sigm_param(:,pixel_off) = [];

% Normalise
ampl_curve_norm = (ampl_curve_sigm-ampl_curve_sigm(1,:))./(ampl_curve_sigm(end,:)-ampl_curve_sigm(1,:));
ampl_curve_norm(:,pixel_off) = [];

figure(1);
subplot(4,4,[11 12 15 16]); hold off;
plot(time/60,ampl_curve_norm,'linewidth',2);
setFigureOptions('Time [min]','Proton concentration [M]',12,'Sigmoidal fitting');

% Plot DNA heat map
figure(1);
subplot(4,4,[2 6]);
surfArrayLW(pixel_dna,'fill',pixel_on);
setFigureOptions('Pixel X','Pixel Y',12,'Max derivative')
colormap jet;

% Plot Ct
figure(1);
subplot(4,4,[3 7]);
surfArrayLW(sigm_param(1,:),'fill',pixel_on(pixel_on_dna));
setFigureOptions('Pixel X','Pixel Y',12,'Max derivative')
colormap jet;

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
for w = 1:nw
    [~,coord_wells_dna{w}] = intersect(pixel_on_dna,coord_wells{w});
    ampl(w) = mean(sigm_param(3,coord_wells_dna{w}))
    Ct(w) = mean(sigm_param(1,coord_wells_dna{w}));
    disp(['Average pH change in well ' num2str(w) ' is ' num2str(ampl(w))])
    disp(['Average Ct in well ' num2str(w) ' is ' num2str(Ct(w)) ' mins'])
end

disp('Over')
end