function y = decimate_avg(x, N, dim)
    mu = floor((N-1) / 2);
    
    y = [];

    for ii = mu+1:N:size(x, dim)-mu
        if dim == 1
            y = [   y;
                    mean(x(ii-mu:ii+mu, :), 1) ];
        elseif dim == 2
            y = [ y, mean(x(:, ii-mu:ii+mu), 2) ];
        end
    end
end