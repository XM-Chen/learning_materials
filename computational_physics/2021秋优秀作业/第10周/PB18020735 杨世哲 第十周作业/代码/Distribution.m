function [x, fx] = Distribution(xi, N)
%xi is an array. The number of its columns is the dimension of the random
%array.
%N is a row vector.
    dim = size(xi, 2);
    x0step = (max(xi) - min(xi))./N;
    xstart = min(xi) - x0step./N;
    xend = max(xi) + x0step./N;
    y = cell(1, dim);
    ShapeArray = cell(1, dim);
    for k = 1:dim
        ShapeArray{k} = 1;
    end
    x = cell(1, dim);
    ind = zeros(size(xi));
    indcff = zeros(size(N, 2), 1);
    
    k = 1;
    y{k} = linspace(xstart(k), xend(k), N(k)+1);
    if dim ~= 1
        ShapeArray{k} = N(k);
        x{k} = reshape((y{k}(2:end) + y{k}(1:end-1))./2, ShapeArray{:});
        ShapeArray{k} = 1;
    else
        x{k} = (y{k}(2:end) + y{k}(1:end-1)).'./2;
    end
    for j = 1:numel(y{k})
        ind(:, k) = (ind(:, k)) + (xi(:, k) > y{k}(j));
    end
    indcff(k) = 1;
    for k = 2:dim
        y{k} = linspace(xstart(k), xend(k), N(k)+1);
        if dim ~= 1
            ShapeArray{k} = N(k);
            x{k} = reshape((y{k}(2:end) + y{k}(1:end-1))./2, ShapeArray{:});
            ShapeArray{k} = 1;
        else
            x{k} = (y{k}(2:end) + y{k}(1:end-1)).'./2;
        end
        for j = 1:numel(y{k})
            ind(:, k) = (ind(:, k)) + (xi(:, k) > y{k}(j));
        end
        ind(:, k) = ind(:, k) - 1;
        indcff(k) = prod(N(1:k-1).');
    end
    
    %Calculate Distribution
    N1 = prod(N.');
    fx = zeros(N1, 1);%Note : high-dimention random array may lead to a very large array.
    ind1 = ind*indcff;
    for k = 1:numel(ind1)
        fx(ind1(k)) = fx(ind1(k)) + 1;
    end
    x1step = (xend - xstart)./N;
    NormFac = sum(fx, 'all').*prod(x1step.');
    fx = fx./NormFac;
    if ~isscalar(N)
        fx = reshape(fx, N);
    end
end