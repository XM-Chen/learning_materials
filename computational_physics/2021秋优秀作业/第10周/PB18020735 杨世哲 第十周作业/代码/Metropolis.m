function x = Metropolis(N, dim, p, New, x0)
%N is the length of the Markov chain.
%dim is the dimension of the random array.
%p is the function handle of the distribution function.
%New is the function for generating x_new.
    x = zeros(N, dim);
    if nargin == 4
        x(1, :) = New(rand(1, dim));
    elseif nargin == 5
        x(1, :) = New(x0);
        xi = rand;
        if xi >= p(x(1, :))/p(x0)
            x(1, :) = x0;
        end
    end
    for k = 2:N
        x(k, :) = New(x(k-1, :));
        xi = rand;
        if xi >= p(x(k, :))/p(x(k-1, :))
            x(k, :) = x(k-1, :);
        end
    end
end