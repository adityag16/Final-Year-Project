function func = denormaliseFct(func_norm,range,func_ref)
%normaliseFct performs linear normalisation of func based on range

func = func_norm*(func_ref(range(2))-func_ref(range(1))) + func_ref(range(1));

end