function setFigureOptions(txt_xlabel,txt_ylabel,val_fontsize,txt_title,val_xlim,val_ylim)

xlabel(txt_xlabel);
ylabel(txt_ylabel);
if nargin > 2
    set(gca,'fontsize',val_fontsize);
end
if nargin > 3
    title(txt_title);
end
if nargin > 4
    xlim(val_xlim);
end
if nargin > 5
    ylim(val_ylim);
end

end

