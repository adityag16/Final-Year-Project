function varargout = initialiseArray(size)
%Initialises all argument to a cell of size N

varargout = cell(nargout,1);
for i = 1:nargout
    varargout{i} = zeros(size);
end

end