function xi = StableMetropolis(N, dim, p, New, Dist, RelTol)
%Generate stable distribution.
%N is the length of the Markov chain.
%p is the function handle of the distribution function.
%New is the function for generating x_new.
%Fitting is a function calculating the parameters of sample distribution.
%Dist is the analytical aim distribution function.
    DistPoints = 80.*ones(1, dim);
    if nargin == 5
        RelTol = 5e-2;
    end
    MaxIterTime = 5;
    Nx = N;
    %-------------First-time generating-------------%
    IterTime = 0;
    xi = Metropolis(Nx, dim, p, New);
    [x, fx] = Distribution(xi, DistPoints);
    fx_check = Dist(x{:});
    fxmax = max(fx_check, [], 'all');
    sig = sum(abs(fx_check - fx) > RelTol.*fxmax, 'all');
    %-------------Iterarion to get stable distribution-------------%
    while sig ~= 0
        IterTime = IterTime + 1;
        if IterTime > MaxIterTime
            warning('Unable to generate a stable distribution');
            break;
        end
        Nx = N.*(2.^IterTime);
        xi = [xi; Metropolis(Nx, dim, p, New, xi(end, :))];
        [x, fx] = Distribution(xi, DistPoints);
        fx_check = Dist(x{:});
        fxmax = max(fx_check, [], 'all');
        sig = sum(abs(fx_check - fx) > RelTol.*fxmax, 'all');
    end
end