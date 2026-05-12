function f_prime = numerical_derivative(f, x, h)
    % numerical_derivative computes the 5-stencil derivative of the scalar function f around x, by step h
    %
    % Arguments:
    %   f: the scalar function f(x)
    %   x: the point around which to compute the derivative
    %   h: the step-size (suggested: 1e-3)
    %
    % See: https://en.wikipedia.org/wiki/Five-point_stencil#1D_first_derivative

    f_prime = (f(x - 2*h) - 8*f(x - h) + 8*f(x + h) - f(x + 2*h)) / (12 * h);
end