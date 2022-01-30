function plotValue(idx,x,y,XY,options)
% Plots a cute value with a y width depending on the data

if nargin<5
    options = '';
end

for i = 1:length(idx)
    if strcmp(XY,'X')
        min_val = y(idx(i))-std(y);
        max_val = y(idx(i))+std(y);
        plot(x(idx(i))*ones(1,100),linspace(min_val,max_val,100),options);
    elseif strcmp(XY,'Y')
        min_val = min(x);
        max_val = max(x);
        plot(linspace(min_val,max_val,100),idx(i)*ones(1,100),options);
    end    
end
end