function [linsys, x_bar, u_bar] = linearize_4tanks(sys, y_bar)
    % linearize_4tanks  Computes the equilibrium of the (non-nominal) Quadruple Tank system and its linearization

    % Make sure that the system contains all the required fields
    so = check_system_params(sys);

    % z = [ x, u ] (Independent variables)

    % Initial guess: all levels to 10cm, inputs to 6V
    z0 = [ 10 * ones(4, 1); 6 * ones(2, 1) ];
    
    % Lower bounds of optimization variables
    z_lb = [ so.x_min * ones(4, 1); 
             so.u_min * ones(2, 1);  ];

    % Upper bounds of optimization variables
    z_ub = [ so.x_max * ones(4, 1); 
             so.u_max * ones(2, 1) ];
    
    % Setup the nonlinear constraints of the opt. problem (i.e., the dynamics)
    nlconstr = @(z) build_nonlinear_constraints(z, so, y_bar);
    
    % The cost function - Ensure that h1 and h3 are as small as possible
    nlcost = @(z) z(1) + z(3);
    
    % Solve the optimization problem and retrieve x_bar and u_bar
    opts = optimoptions(@fmincon,'Algorithm','sqp');
    z_opt = fmincon(nlcost, z0, [], [], [], [], z_lb, z_ub, nlconstr, opts);
    x_bar = reshape(z_opt(1:4), [4, 1]);
    u_bar = reshape(z_opt(5:6), [2, 1]);

    % Linearize the system
    % d f(x, u) / dx
    a11 = - so.A(1) * sqrt(so.g) / (so.S * sqrt(2*x_bar(1)));
    a21 = so.A(1) * sqrt(so.g) / (so.S * sqrt(2*x_bar(1)));
    a22 = - so.A(2) * sqrt(so.g) / (so.S * sqrt(2*x_bar(2)));
    a33 = - so.A(3) * sqrt(so.g) / (so.S * sqrt(2*x_bar(3)));
    a43 = so.A(3) * sqrt(so.g) / (so.S * sqrt(2*x_bar(3)));
    a44 = - so.A(4) * sqrt(so.g) / (so.S * sqrt(2*x_bar(4)));
    
    % d f(x, u) / du (Computed numerically using the 5-point stencil method)
    b31 = numerical_derivative(@(v) so.k1(v) * so.gamma1(v) * v, u_bar(1), 1e-3) / so.S;
    b12 = numerical_derivative(@(v) so.k2(v) * so.gamma2(v) * v, u_bar(2), 1e-3) / so.S;
    b21 = numerical_derivative(@(v) so.k1(v) * (1 - so.gamma1(v)) * v, u_bar(1), 1e-3) / so.S;
    b42 = numerical_derivative(@(v) so.k2(v) * (1 - so.gamma2(v)) * v, u_bar(2), 1e-3) / so.S;

    A = [   a11,    0,      0,      0;
            a21,    a22,    0,      0;
            0,      0,      a33,    0;
            0,      0,      a43,    a44 ];

    B = [   0,      b12;
            b21,    0;
            b31,    0;
            0,      b42 ];

    C = [   0,      1,      0,      0;
            0,      0,      0,      1; ];

    D = zeros(2, 2);

    linsys = ss(A, B, C, D);
end

function [ nlcnstr_le, nlcnstr_eq ] = build_nonlinear_constraints(z, sys, y_bar)
    % build_nonlinear_constraints  Constructs the system's nonlinear constraints in a fmincon-friendly way
    % Set of constraintz h(z) <= 0
    nlcnstr_le = [];

    % Set of constraints h(z) == 0 (DYNAMICS)
    nlcnstr_eq = [ dynamics_f(z(1:4), z(5:6), sys);
                   y_bar - dynamics_g(z(1:4)) ];
end


function x_dot = dynamics_f(x, u, sys)
    % dynamics_f Implements the system's dynamics

    x_dot = zeros(4, 1);

    x_dot(1) = - sys.A(1) / sys.S * sqrt(2 * sys.g * x(1)) + (sys.k2(u(2)) * sys.gamma2(u(2)) / sys.S) * u(2);
    x_dot(2) = + sys.A(1) / sys.S * sqrt(2 * sys.g * x(1)) - sys.A(2) / sys.S * sqrt(2 * sys.g * x(2)) ...
                + (sys.k1(u(1)) * (1 - sys.gamma1(u(1))) / sys.S) * u(1);
    x_dot(3) = - sys.A(3) / sys.S * sqrt(2 * sys.g * x(3)) + (sys.k1(u(1)) * sys.gamma1(u(1)) / sys.S) * u(1);
    x_dot(4) = + sys.A(3) / sys.S * sqrt(2 * sys.g * x(3)) - sys.A(4) / sys.S * sqrt(2 * sys.g * x(4)) ...
                + (sys.k2(u(2)) * (1 - sys.gamma2(u(2))) / sys.S) * u(2);
end

function y = dynamics_g(x)
    % dynamics_g Implements the system's output transformation
    y = [ x(2); x(4) ];
end

function sys_ok = check_system_params(sys)
    % check_system_params Fills possible missing paramers with the nominal values
  
    sys_ok = struct;
    sys_ok.S = getfield_default(sys, 'S', (4.44 / 2)^2 * pi);
    sys_ok.A = getfield_default(sys, 'A', [0.0792, 0.2425, 0.0792, 0.2425]);
    sys_ok.g = getfield_default(sys, 'g', 981.0);
    sys_ok.k1 = getfield_default(sys, 'k1_fcn',  @(V) 3.3);
    sys_ok.k2 = getfield_default(sys, 'k2_fcn',  @(V) 3.3);
    sys_ok.gamma1 = getfield_default(sys, 'gamma1_fcn',  @(V) 0.36);
    sys_ok.gamma2 = getfield_default(sys, 'gamma2_fcn',  @(V) 0.36);
    sys_ok.u_min = getfield_default(sys, 'u_min', 0.0);
    sys_ok.u_max = getfield_default(sys, 'u_max', 12.0);
    sys_ok.x_min = getfield_default(sys, 'x_min', 0.1);
    sys_ok.x_max = getfield_default(sys, 'x_max', 30.0);
end

function val = getfield_default(s, field, default)
    if isfield(s, field)
        val = s.(field);
    else
        val = default;
    end
end