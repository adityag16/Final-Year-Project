 function varargout = cellparam(options,str,def)
%cellparam finds the parameters following a the presence of strings str in
% cell options.

% Initialisation
N = length(str);
param = zeros(N,1);
if nargin<3
    def = zeros(N,1);
end
varargout = cell(N,1);

% Execution
for i = 1:N
    if cellContainsStr(options,str{i})
        [~,idx] = cellContainsStr(options,str{i});
        varargout{i} = options{idx+1};
    else
        varargout{i} = def(i);
    end
end