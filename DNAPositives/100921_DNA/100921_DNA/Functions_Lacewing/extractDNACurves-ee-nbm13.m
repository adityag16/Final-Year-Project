function [time,out_smooth,ampl_curve_V,ampl_curve_pH,ampl_curve_lin,ampl_curve_sigm,...
          exp_comp,sigm_param,pixel_on,pixel_on_dna,xm,xi,xf,Ct,ampl] ...
          = extractDNACurves(file_data,file_data_all,start_exp,end_exp,pH_init,options)
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
%       - 'spreadpixels' : picks a sub-array of rxc pixels within the array to
%         debug output. [r c] must be the next argument.
%       - 'custompixels' : considers a preloaded vector of pixel 
%         coordinates pixel_coord to apply the function to. pixel_coord
%         must be the next argument.
%       - 'skipaverage' : do not run algo on average output.
%       - 'showallpixels' : show curves for each pixel.
%       - 'onlyaverage' : only run algo on average.
%       - 'wells' : if experiment has two wells.
      
% Constants
[ROWS,COLS,np,coord_chem,ROWS_temp,COLS_temp,npt,coord_temp,coord_peltier] = getConstantsLW;

% Algorithm parameters
sens = 12;
minute = 60;
smooth_factor = 40;

% Option
if nargin<6
    options = {};
end
if cellContainsStr(options,'spreadpixels')
    [~,idx] = cellContainsStr(options,'spreadpixels');
    pixel_spread = options{idx+1};
end
if cellContainsStr(options,'custompixels')
    [~,idx] = cellContainsStr(options,'custompixels');
    pixel_custom = options{idx+1};
end
if cellContainsStr(options,'wells')
    [~,idx] = cellContainsStr(options,'wells');
    nw = options{idx+1};% Number of wells
else
    nw=1;
end

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
if cellContainsStr(options,'startchem')
    en_temp = 0;
else
    en_temp = 1;
end

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
[pixel_on, pixel_off, npon,npoff] = removeOutPixelsLW(out,pixel_on_Vref,x_stable);

% Read a set of r x c pixels throughout the array
if cellContainsStr(options,'spreadpixels')
    pixel_on = generateSubArrayLW(pixel_spread);
    [~,idx] = intersect(pixel_on,coord_temp);
    if ~isempty(idx)
        pixel_on(idx) = pixel_on(idx)+1;
    end
    npon = length(pixel_on);
    pixel_off = setdiff(1:np,pixel_on);
    npoff = length(pixel_off);
elseif cellContainsStr(options,'custompixels')
    pixel_on = pixel_custom;
    npon = length(pixel_on);
    pixel_off = setdiff(1:np,pixel_on);
    npoff = length(pixel_off);
end

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

% Average of array
if nw==1
    out_smooth_av = mean(out_smooth,2);
else
    out_smooth_av = zeros(length(time),nw)
    thr = 40;
    [coord_wells,wells] = divideArrayLW(thr);
    for w = 1:nw
        [~,coord_wells{w}] = intersect(pixel_on,coord_wells{w});

        out_smooth_av(:,w) = mean(out_smooth(:,coord_wells{w}),2);%mean(out_smooth(:,coord_wells{1}),2);
    end
end

% Plot average of all pixels
figclf(1);
ax1 = subplot(4,4,1);
plot(time/60,out_smooth_av);
title('All pixels');
setFigureOptions('','Voltage [mV]',12,'All pixels')
ax2 = subplot(4,4,5);
plot(time/60,smooth(Vtemp,smooth_factor),'linewidth',2);
linkaxes([ax1,ax2],'x');
setFigureOptions('Time [min]','Voltage [mV]',12)

if cellContainsStr(options,'subtract')
    nw = 1;
    out_smooth_av = mean(out_smooth(:,coord_wells{2}),2)-mean(out_smooth(:,coord_wells{1}),2);
end
    
for w=1:nw
    % Find maximum of derivative
    x0(w) = findDNAPeakLW(out_smooth_av(:,w),time,smooth_factor,x_stable,[time(x_stable) time(end)-3*60]);
    disp(['Mean derivative maximum at ' num2str(time(x0(w))/60) ' min']);
end

% Cut experiment to double xm
% out_smooth_av = out_smooth_av(1:min(length(time),2*x0));
% time = time(1:min(length(time),2*x0));

% Plot
figure(1); 
subplot(4,4,[2 6]); hold all;
plot(time/60,out_smooth_av,'linewidth',2);
setFigureOptions('Time [min]','Voltage [mV]',12,'Average output')

if ~(cellContainsStr(options,'spreadpixels')||cellContainsStr(options,'skipaverage'))
    for w = 1:nw
        % Plot derivative max
        figure(1); 
        subplot(4,4,[2 6]);
        plotValue(x0(w),time/60,out_smooth_av(:,w),'X');
        
        % Find amplification region for average
        [xi_av(w),xf_av(w),out_exp_i(:,w),out_exp_f(:,w)] = findDNABoundsLW(out_smooth_av(:,w),time,x_stable,x0(w),smooth_factor);
        % Plot
        figure(1); subplot(4,4,[2 6]); hold all;
        plotValue([x0(w),xi_av(w),xf_av(w)],time/60,out_smooth_av(:,w),'X');
%         figure(1); subplot(4,4,[3 7]); hold all;
%         plotValue([x0(w),xi_av(w),xf_av(w)],time/60,out_smooth_av(:,w),'X');
    %     plot(time/60,out_exp_i,'--','linewidth',2);
    %     plot(time/60,out_exp_f,'--','linewidth',2);

        % Drift compensation
        [ampl_curve_V_av(:,w),exp_comp_av(:,w)] = compensateDriftLW(time,out_smooth_av(:,w),xi_av(w),xf_av(w),out_exp_i(:,w),out_exp_f(:,w),x_stable);
        figure(1); 
        subplot(4,4,[2 6]);
        plot(time/60,exp_comp_av(:,w),'--','linewidth',2);
        subplot(4,4,3); hold all;
        plot(time/60,ampl_curve_V_av(:,w),'linewidth',2);
        plotValue([x0(w),xi_av(w),xf_av(w)],time/60,ampl_curve_V_av(:,w),'X');
        setFigureOptions('Time [min]','Voltage [mV]',12,'Voltage curve')

        % Linearise amplification curvess
        [ampl_curve_lin_av(:,w),ampl_curve_pH_av(:,w)] = linAmplCurveLW(ampl_curve_V_av(:,w),x_stable,pH_init,smooth_factor);

        subplot(4,4,7); hold all;
        plot(time/60,ampl_curve_lin_av(:,w),'linewidth',2);
        setFigureOptions('Time [min]','Voltage [mV]',12,'Linear curve')

        % Sigmoidal fitting
        [ampl_curve_sigm_av(:,w) sigm_param_av] = fitSigmoidLW(ampl_curve_lin_av(:,w),time,x_stable);
        plot(time/60,ampl_curve_sigm_av(:,w),'--','linewidth',2);
        
        subplot(4,4,[4 8]);
        surfArrayLW(ones(npon,1),'fill',pixel_on);
        setFigureOptions('Pixel X','Pixel Y',12,'On pixels');
        
        thr = ampl_curve_sigm_av(1,w) + .2*(ampl_curve_sigm_av(end,w)-ampl_curve_sigm_av(1,w));
        [~,idx] = min(abs(ampl_curve_sigm_av(:,w)-thr));
        TTP = time(idx)/60
    end
elseif cellContainsStr(options,'spreadpixels')
    xi_av = 1;
    xf_av = length(time);
    % If spread pixels, display the chosen pixels
    subplot(4,4,[4 8]);
    surfArrayLW(ones(npon,1),'fill',pixel_on);
    setFigureOptions('Pixel X','Pixel Y',12,'Sub-array of pixels');
    colormap parula;
end

if cellContainsStr(options,'onlyaverage')
    out_smooth=out_smooth_av;
    ampl_curve_V=ampl_curve_V_av;
    ampl_curve_pH=ampl_curve_pH_av;
    ampl_curve_lin=ampl_curve_lin_av;
    ampl_curve_sigm=ampl_curve_sigm_av;
    exp_comp=exp_comp_av;
    sigm_param=sigm_param_av;
    pixel_on=0;
    pixel_on_dna=0;
    xm=x0;
    xi=xi_av;
    xf=xf_av;
    Ct = TTP;
    ampl = 0;
    return;
end

%% Then for each pixel

% Constants
npon = length(pixel_on_dna);

% Find max derivatives
range = [time(x_stable) time(end)-3*60];
% range = [time(xi_av)-3*60,time(xf_av)+3*60];
[xm,tm] = findDNAPeakLW(out_smooth,time,4*smooth_factor,x_stable,range);

% Sort positives and negatives - 0 unstable, 1 positive, 2 negative
pixel_dna = ones(size(pixel_on));
pixel_dna(union(find(tm<range(1)+2*60),find(tm>range(2)-2*60))) = 2;

% Plot all max derivatives
figure(1);
subplot(4,4,[4 8]);
surfArrayLW(tm/60,'fill',pixel_on_dna);
setFigureOptions('Pixel X','Pixel Y',12,'Max derivative')
colormap jet;

% Plot DNA heat map
figure(1);
subplot(4,4,[9 13]);
surfArrayLW(pixel_dna,'fill',pixel_on);
setFigureOptions('Pixel X','Pixel Y',12,'Max derivative')
colormap jet;

% Initialisation
[out_exp_i,out_exp_f,ampl_curve_V,ampl_curve_pH,ampl_curve_lin,exp_comp,ampl_curve_sigm] = initialiseArray([n,npon]);
[xi,xf,error_fit] = initialiseArray([npon,1]);
sigm_param= zeros(4,npon);
% Waitbar
f = waitbar(0,'Processing amplification curves');
% p = find(pixel_on == 78*36+10);
for p = 1:npon
    % Suspicion of positive
    if pixel_dna(p) == 1
        % Find amplification region for average
        [xi(p),xf(p),out_exp_i(:,p),out_exp_f(:,p)] = ...
            findDNABoundsLW(out_smooth(:,p),time,x_stable,xm(p),smooth_factor);

%         figclf(3); hold all;
%         plot(time/60,out_smooth(:,p));
%         plotValue([xm(p) xi(p) xf(p)],time/60,out_smooth(:,p),'X');
%         plot(time/60,out_exp_i(:,p));
%         plot(time/60,out_exp_f(:,p));
        
        % Drift compensation
        [ampl_curve_V(:,p),exp_comp(:,p)] = ...
            compensateDriftLW(time,out_smooth(:,p),xi(p),xf(p),...
            out_exp_i(:,p),out_exp_f(:,p),x_stable);

        % Linearisation
        [ampl_curve_lin(:,p),ampl_curve_pH(:,p)] = linAmplCurveLW(ampl_curve_V(:,p),x_stable,pH_init,smooth_factor);

        % Sigmoidal fitting
        [ampl_curve_sigm(:,p),sigmoid_fit,sigmoid_gof,sigmoid_scale] = fitSigmoidLW(ampl_curve_lin(:,p),time,x_stable);
        error_fit(p) = sigmoid_gof.rsquare;

        % Sigmoidal feature extraction
        if error_fit(p)>.95 && (ampl_curve_lin(end,p)>ampl_curve_lin(1,p))
            sigm_param(:,p) = extractSigmParamLW(sigmoid_fit,sigmoid_scale,time,pH_init);
        else
            pixel_dna(p) = 2; % might be negative
        end
    end
    
    if pixel_dna(p) == 2
        % Check if negative by fitting
        [~,~,rsquare_neg] = compensateExp2LW(out_smooth(:,p),time,[time(x_stable),time(end)],x_stable);
        if rsquare_neg<.95 % unstable output
            pixel_dna(p) = 0;
        end
        ampl_curve_V(:,p) = 500*ones(size(out_smooth(:,p)));
        ampl_curve_pH(:,p) = pH_init*ones(size(out_smooth(:,p)));
        ampl_curve_lin(:,p) = 10^(-pH_init)*ones(size(out_smooth(:,p)));
        ampl_curve_sigm(:,p) = ampl_curve_lin(:,p);
    end
    
    if cellContainsStr(options,'spreadpixels') || cellContainsStr(options,'showallpixels')
        if mod(p,1)==0
            % Highlight pixel being read
            subplot(4,4,[4 8]);
            vect = .5*ones(npon,1); vect(p)=1;
            surfArrayLW(vect,'fill',pixel_on);
            colormap parula;

            % Display output
            subplot(4,4,[9 10]);
            plot(time/60,out_smooth(:,p)); hold on;
            plot(time/60,out_exp_i(:,p));
            plot(time/60,out_exp_f(:,p));
%             plotValue([xm(p),xi(p),xf(p)],time/60,out_smooth(:,p),'X');
            setFigureOptions('','Voltage [mV]',12,'Raw output')
            hold off;

            subplot(4,4,[13 14]);
            plot(time/60,ampl_curve_V(:,p)); hold on;
%             plotValue([xm(p),xi(p),xf(p)],time/60,ampl_curve_V(:,p),'X');
            setFigureOptions('Time [min]','Voltage [mV]',12,'Voltage curve')
            hold off;

            subplot(4,4,[11 12]);
            plot(time/60,ampl_curve_lin(:,p));
            setFigureOptions('','Voltage [mV]',12,'Linear curve')

            subplot(4,4,[15 16]);
            plot(time/60,ampl_curve_sigm(:,p));
            setFigureOptions('Time [min]','Proton concentration',12,'Sigmoid')
        
            waitforbuttonpress;
        end
    end
       
    waitbar(p/npon, f,sprintf('Pixel number %f / %f',p,npon));
end
close(f);

% Post-filtering
pixel_dna(sigm_param(1,:)<4) = 2;
pixel_dna(sigm_param(3,:)>1) = 2;
pixel_dna(ampl_curve_sigm(end,:)-ampl_curve_sigm(1,:)>1e-7) = 2;

% Normalise
ampl_curve_norm = (ampl_curve_sigm-ampl_curve_sigm(1,:))./(ampl_curve_sigm(end,:)-ampl_curve_sigm(1,:));

figure(1);
subplot(4,4,[10 11 14 15]); hold off;
plot(time/60,ampl_curve_norm(:,(pixel_dna==1)),'linewidth',2);
setFigureOptions('Time [min]','Proton concentration [M]',12,'Sigmoidal fitting');

% Plot DNA heat map
figure(1);
subplot(4,4,[9 13]);
surfArrayLW(pixel_dna,'fill',pixel_on);
setFigureOptions('Pixel X','Pixel Y',12,'Max derivative')
colormap jet;

% Plot Ct
figure(1);
subplot(4,4,[12 16]);
surfArrayLW(sigm_param(1,pixel_dna==1),'fill',pixel_on(pixel_dna==1));
setFigureOptions('Pixel X','Pixel Y',12,'Max derivative')
colormap jet;



% Feature extraction
colormap jet;
feature_label = {'Ct','Slope','Amplitude','Asymmetry'};
figclf(2);
for i = 1:4
    subplot(2,2,i);
    surfArrayLW(sigm_param(i,pixel_dna==1),'fill',pixel_on(pixel_dna==1));
    setFigureOptions('Pixel X','Pixel Y',12,feature_label{i})
    view(2);
    colorbar;
end

if cellContainsStr(options,'wells')
    [~,idx] = cellContainsStr(options,'wells');
    nw = options{idx+1};% Number of wells
else
    nw=1;
end

% Average amplitude and Ct
pixel_on_dna = find(pixel_dna==1);
Ct = zeros(nw,1);
ampl = zeros(nw,1);
for w = 1:nw
    [~,coord_wells_dna{w}] = intersect(pixel_on(pixel_on_dna),pixel_on(coord_wells{w}));
    ampl(w) = mean(sigm_param(3,pixel_on_dna(coord_wells_dna{w})));
    Ct(w) = mean(sigm_param(1,pixel_on_dna(coord_wells_dna{w})));
    disp(['Average pH change in well ' num2str(w) ' is ' num2str(ampl(w))])
    disp(['Average Ct in well ' num2str(w) ' is ' num2str(Ct(w)) ' mins'])
end

disp('Over')
end