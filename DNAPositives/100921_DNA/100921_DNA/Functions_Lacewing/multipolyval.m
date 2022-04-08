function F = multipolyval(P, X)

    n_signals = size(P, 2);

    if isrow(X)
       X = repmat(X', 1, n_signals);
    elseif iscolumn(X)
       X = repmat(X, 1, n_signals);
    else
        assert(size(X, 2) == n_signals, 'Something is wrong...');
    end

    F = zeros(size(X, 1), n_signals);

    for i = 1:n_signals
        F(:, i) = polyval(P(:, i), X(:, i));
    end

end

