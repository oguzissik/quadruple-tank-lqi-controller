function y = fill_nans(x, dim)
    y = x;
    maxels = max(x, [], dim);
    
    if dim == 1
        maxels = repmat(maxels, size(x, 1), 1);
    else
        maxels = repmat(maxels, size(x, 2), 2);
    end

    y(isnan(y)) = maxels(isnan(y));

end

    
