function func_norm = normaliseFct(func,range)
%normaliseFct performs linear normalisation of func based on range

func_norm = (func - func(range(1)))./...
    (func(range(2))-func(range(1)));

end