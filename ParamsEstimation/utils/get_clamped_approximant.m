function fcn = get_clamped_approximant(x, y, x_int, y_int, type)
    [xData, yData] = prepareCurveData(x, y);
    if strcmp(type, 'spline')
        ft = fittype( 'smoothingspline' );
        opts = fitoptions( 'Method', 'SmoothingSpline' );
        opts.Normalize = 'on';
        opts.SmoothingParam = 0.9425;
        [fcn_ucp, ~] = fit(xData, yData, ft, opts);
    elseif strcmp(type, 'logistic')
        [~, idx] = min(x);
        y_min = y(idx);

        equation = sprintf('a*logistic(b*(x+c))-a*logistic(b*(%d+c))+%d', x_int(1), y_min);
        ft = fittype(equation, 'independent', 'x', 'dependent', 'y', 'coefficients', {'a', 'b', 'c'});
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [0.75, 0.45, 0.0];
        
        % Fit model to data.
        [fcn_ucp, ~] = fit( xData, yData, ft, opts );
    elseif strcmp(type, 'sin')
        ft = fittype('sin2');
        opts = fitoptions('Method', 'NonlinearLeastSquares');
        opts.Display = 'Off';
        opts.Robust = 'LAR';
        opts.Lower = [-Inf 0 -Inf -Inf 0 -Inf];
        opts.StartPoint = [5.34, 0.29, -1.2, 2.24, 0.57, -0.71];

        % Fit model to data.
        [fcn_ucp, ~] = fit(xData, yData, ft, opts);
    else
        error('Unsupported!')
    end

    fcn = @(x) clamped_function(fcn_ucp, x, x_int(1), x_int(2), y_int(1), y_int(2));
end

function y = clamped_function(fun, x, x_min, x_max, y_min, y_max)
    if isempty(y_min) || isnan(y_min)
        y_min = fun(x_min);
    end
    if isempty(y_max) || isnan(y_max)
        y_max = fun(x_max);
    end

    y = y_min .* (x < x_min) + y_max .* (x > x_max) ...
            + reshape(fun(x(:)), size(x, 1), size(x, 2)) .* (1-(x < x_min)) .* (1-(x > x_max));
end