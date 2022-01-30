function P = multipolyfit(X, Y, order)

    assert(all(size(X)==size(Y)), 'Get it right son!')

    [n_observations, n_signals] = size(X);

    P = zeros(order+1, n_signals);

    for i = 1:n_signals
        P(:, i) = polyfit(X(:, i), Y(:, i), order);
    end

end

