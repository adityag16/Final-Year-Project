function [bool,idx] = cellContainsStr(C,str)
%cellContainsStr checks if cell contains the string str

    % Look for string
    idx = cellfun(@(s) ~isempty(strfind(str, s)), C);
    idx = find(idx==1);
    bool = ~isempty(idx);

end

