function Ct = findCt(curve,time,thr)

idx_Ct = find(curve<thr);
idx_Ct = idx_Ct(end);
if (idx_Ct > 1) && (idx_Ct < length(time))
    slope = diff(curve(idx_Ct:idx_Ct+1))/diff(time(idx_Ct:idx_Ct+1));
    Ct = time(idx_Ct) + ...
        (thr-curve(idx_Ct))/(slope);
else
    Ct = 0;
end

end