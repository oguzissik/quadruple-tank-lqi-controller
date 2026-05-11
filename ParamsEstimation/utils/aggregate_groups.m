function y = aggregate_groups(groups, x, varargin)
    aggregate_dim = 2;

    if numel(varargin) >= 1
        aggregate_dim = varargin{1};
    end

    n_groups = numel(unique(groups));
    
    if aggregate_dim == 1
        y = zeros(n_groups, size(x, 2));
        for ii=1:size(x, 2)
            y(:, ii) = accumarray(groups(:), x(:, ii), [], @mean);
        end
    else
        y = zeros(size(x, 1), n_groups);
        for ii=1:size(x, 1)
            y(ii, :) = accumarray(groups(:), x(ii, :).', [], @mean).';
        end
    end

end