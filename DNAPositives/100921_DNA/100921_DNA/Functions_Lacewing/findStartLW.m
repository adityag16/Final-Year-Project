function [range,x_stable,t_stable] = findStartLW(time,out,en_temp,start_exp,end_exp,nw)
%findStartLW identifies the time of injection and the time of settling to
% start performing data processing.q

% Constants
[ROWS,COLS,np,coord_chem,ROWS_temp,COLS_temp,npt,coord_temp] = getConstantsLW;

% Finding the range
range = [0 0];
% End of measurement
range(2) = length(time);
% Approximate time of loading
% Variable definition
if en_temp
    out_temp = mean(out(coord_temp,:),1);
else
    out_temp = mean(out(coord_chem,:),1);
end
out_temp_diff = diff((out_temp));
out_temp_diff_smooth = diff(smooth(out_temp));

% % Figure
% figure(1);
% ax1 = subplot(2,1,1);
% plot(time,out_temp);
% ax2 = subplot(2,1,2);
% plot(time(1:end-1),out_temp_diff);
% linkaxes([ax1,ax2],'x');

% Define bounds for the start of the experiment
x_start = time_to_idx(time, start_exp);
x_start = x_start(1):x_start(2);

if nw==1
    if start_exp(1)==start_exp(2)
        range(1) = time_to_idx(time,start_exp(1));
    else
        % Define bounds for the start of the experiment
        x_start = time_to_idx(time, start_exp);
        x_start = x_start(1):x_start(2);
        % Find the time at which temperature drops as the latest after 
        idx        = find(out_temp_diff(x_start)<-.7) + x_start(1)-1;
        if en_temp==0
            idx        = find(out_temp_diff(x_start)<-.1) + x_start(1)-1;
        end
        % If empty, loosen constraints
        if isempty(idx)
            idx        = find(out_temp_diff(x_start)<-.3) + x_start(1)-1;
        end
        if length(idx)>1
            idx_smooth = find(out_temp_diff_smooth(x_start)<-.35) + x_start(1)-1;
            % If empty, loosen constraints
            if isempty(idx_smooth)
                idx_smooth = find(out_temp_diff_smooth(x_start)<-.1) + x_start(1)-1;
            end
            [~,idx_min] = min(out_temp_diff(idx_smooth));
            id = find(idx>idx_smooth(idx_min)-3);
            idx = idx(id(1));
        end
        %diff_idx = diff(idx)
        %idx = idx(end);
        % If there are many points, take the closest to 700s

        % Find the first point where the temperature increases again
        idx2 = find(diff(out_temp(idx:x_start(end)))>0);
        if ~isempty(idx2)
            range(1) = idx + idx2(1)-1+2;
        else
            range(1) = idx+2;
        end
        
        out_temp(range(1):end) = smooth(out_temp(range(1):end));
    end
    disp(['Detected loading time ' num2str(time(range(1))) ' s'])
    
    % Calculate derivative of out_temp
    out_temp_diff = abs(smooth(diff(out_temp),50));
    figure(6); 
    subplot(2,1,1);
    plot(time,out_temp);
    subplot(2,1,2);
    plot(time(1:end-1),out_temp_diff);
    xlim([time(range(1)) time(end-1)]);
    [out_temp_diff_max,idx_max] = max(out_temp_diff(range(1):end));
    idx_settled = find(out_temp_diff(range(1):end)<0.05*out_temp_diff_max);
    idx = find(idx_settled>idx_max);
    idx_settled = idx_settled(idx(1)) + range(1) - 1; % Index of settled value
        
    % Detecting t_stable i.e. time at which the output becomes stable
    % and output curves will be considered as real data
    % Considered as reaching 90% of temperature after 2 mins
    temp_heating(1) = out_temp(range(1));
%     [~,idx] = min(abs(time-(time(range(1))+240)));
    temp_heating(2) = mean(out_temp(idx_settled:idx_settled+20));
    x_stable = find(out_temp(range(1)+idx_max:range(2))>(temp_heating(1)+.9*(temp_heating(2)-temp_heating(1))));
    x_stable = x_stable(1) + range(1) + idx_max - 1;
    t_stable = time(x_stable);
    disp(['Detected settling time ' num2str(t_stable) ' s'])

    range(2) = time_to_idx(time,end_exp);
    x_stable = x_stable - range(1) + 1;
    if en_temp==0 && x_stable < range(1)+10
        x_stable = x_stable+30;
    end
    t_stable = time(x_stable);
elseif nw==2
    % Find coordinates of both halves of the array
    coord_wells = divideArrayLW(40);
    coord_well_temp = cell(nw,1);
    figclf(5);
    for w = 1:nw
        coord_well_temp{w} = intersect(coord_wells{w},coord_temp);
        coord_well_chem{w} = intersect(coord_wells{w},coord_chem);
        out_temp_well(:,w) = mean(out(coord_well_temp{w},:),1);
        out_chem_well(:,w) = mean(out(coord_well_chem{w},:),1);
        diff_out_temp_well(:,w) = diff(mean(out(coord_well_temp{w},:),1));
        diff_out_chem_well(:,w) = diff(mean(out(coord_well_chem{w},:),1));
        
        % Plot
        ax(1)=subplot(4,2,w);
        plot(time,out_chem_well(:,w));
        ax(2)=subplot(4,2,w+2);
        plot(time(2:end),diff_out_chem_well(:,w));
        ax(3)=subplot(4,2,w+4);
        plot(time(3:end),diff(diff_out_chem_well(:,w)));
        ax(3)=subplot(4,2,w+4);
        plot(time,out_temp_well(:,w));
        ax(4)=subplot(4,2,w+6);
        plot(time(2:end),diff_out_temp_well(:,w));
        linkaxes(ax,'x')
        
        [y,x] = findpeaks(abs(diff_out_chem_well(x_start,w)));
        idx = x_start(1) + x(find(y==max(y)));
        
%         idx2 = find(abs(diff_out_chem_well(idx:x_start(2),w))<.1);
%         idx2 = idx(1);
        
        % Find the first point where the temperature increases again
%         thr = .1; % Manual threshold
        thr = mean(abs(diff(diff_out_chem_well(idx+60:min(idx+80,length(time)),w)))); % Automatic threshold
        idx2 = find(abs(diff(diff_out_chem_well(idx:x_start(end),w)))<1.5*thr);
        if ~isempty(idx2)
            range(1) = idx + idx2(1) + 1;
        else
            range(1) = idx;
        end
                
        disp(['Detected loading time for well ' num2str(w) ' is ' ...
            num2str(time(range(1))) ' s'])
    end
    
    % Detecting t_stable i.e. time at which the output becomes stable
    % and output curves will be considered as real data
    % Considered as reaching 90% of temperature after 2 mins
    temp_heating(1) = out_temp_well(range(1),w);
    [~,idx] = min(abs(time-(time(range(1))+200)));
    temp_heating(2) = mean(out_temp_well(idx-5:idx+5,w));
    x_stable = find(out_temp_well(range(1):range(2),w)>(temp_heating(1)+.95*(temp_heating(2)-temp_heating(1))));
    x_stable = x_stable(1) + range(1) - 1;
    t_stable = time(x_stable);
    disp(['Detected settling time ' num2str(t_stable) ' s'])

    range(2) = time_to_idx(time,end_exp);
    x_stable = x_stable - range(1) + 1;
    t_stable = time(x_stable);
   
end


