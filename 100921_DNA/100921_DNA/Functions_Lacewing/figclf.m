function h = figclf(n,opt)
%Creates a figure and clears it
% if opt=1 then hold

if nargin<2
    opt = 0;
end

h = figure(n); clf(n);
if opt
    hold all;
end

end

