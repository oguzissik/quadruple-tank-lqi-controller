addpath ../ParamsEstimation/utils;
load('../ParamsEstimation/pump_params_april.mat');
load('../ParamsEstimation/valve_opening.mat');

S = (4.44 / 2)^2 * pi;
y_target = [10; 10];

sys = struct('A', AA, 'S', S, 'g', 981, 'k1_fcn', k1_fcn, 'k2_fcn', k2_fcn, ...
             'gamma1_fcn', gamma1_fcn, 'gamma2_fcn', gamma2_fcn, ...
             'u_min', 4, 'u_max', 12, 'x_min', 0.1, 'x_max', 30);

% Initial state (for simulation)
h0 = [ 2, 2, 2, 2 ].';

[ linsys, x_bar, u_bar ] = linearize_4tanks(sys, y_target);
h_0 = x_bar;

Gs = tf(linsys);
R1 = pidtune(Gs(1, 1), 'PI', 0.05);
R2 = pidtune(Gs(2, 2), 'PI', 0.05);

% Approximate the maps -> They are implemented via Lookup tables
xv = 0:0.05:18;
yk1 = k1_fcn(xv);
yk2 = k2_fcn(xv);
ygamma1 = gamma1_fcn(xv);
ygamma2 = gamma2_fcn(xv);

% Backwards-decoupler
T12 = minreal(- Gs(1, 2) / Gs(1, 1));
T21 = minreal(- Gs(2, 1) / Gs(2, 2));
TT = [  0,  T12;
       T21, 0 ];

Delta = inv((eye(2) - TT));

L_decentralized = Gs * blkdiag(R1, R2);
L_decoupled = Gs * Delta * blkdiag(R1, R2);

% uncomment below to see singular values plot
% figure;
% sigmaplot(L_decentralized, L_decoupled, Gs);
% legend('Decentralized', 'Decoupled', 'Open Loop');
% set(gca, 'FontSize', 16); version error , de-comment in 2021b

%% Discretize Decentralized/Decoupled

% R1d = c2d(R1, 0.01, 'tustin');
% R2d = c2d(R2, 0.01, 'tustin');
% T12d = c2d(T12, 0.01, 'tustin');
% T21d = c2d(T21, 0.01, 'tustin');
% save('DiscretizedController', 'R1d', 'R2d', 'T12d', 'T21d');


%% LQI Design

% Augmented system: add integral states on outputs (h2 and h4)
% State vector: [h1_tilde, h2_tilde, h3_tilde, h4_tilde, e1, e2]
% e1 = integral(h2_ref - h2), e2 = integral(h4_ref - h4)
A_aug = [linsys.A, zeros(4, 2);
         -linsys.C, zeros(2, 2)];

B_aug = [linsys.B;
         zeros(2, 2)];

% Q matrix: penalize state deviations and integral errors
% Higher weight on h2 and h4 (outputs we care about)
% Higher weight on integral states to ensure zero steady-state error
Q = diag([0.1, 1.2, 0.1, 1.2, 0.025, 0.025]);

% R matrix: penalize control effort
% Larger R = more conservative, smaller voltages
R_lq = diag([2.9, 2.9]);

% Solve LQR - find optimal gain matrix K
K_aug = lqr(A_aug, B_aug, Q, R_lq);

% Split K into state feedback and integral feedback parts
Kx = K_aug(:, 1:4);   % 2x4 - state feedback gain
Ke = K_aug(:, 5:6);   % 2x2 - integral feedback gain

disp('LQI gains computed:');
disp('Kx = '); disp(Kx);
disp('Ke = '); disp(Ke);


%% Observer Design

% Observer poles should be 3-5x faster than controller poles
% Controller poles - eigenvalues of (A - B*Kx)
ctrl_poles = eig(linsys.A - linsys.B * Kx);
disp('Controller poles:'); disp(ctrl_poles);

% Place observer poles 5x faster (more negative)
obs_poles = 5 * ctrl_poles;
disp('Observer poles:'); disp(obs_poles);

% Compute observer gain using pole placement
L_obs = place(linsys.A', linsys.C', obs_poles)';
disp('Observer gain L:'); disp(L_obs);

%% Discretize LQI
Ts = 0.01;

% Discretize augmented system (for LQI gains)
sys_aug_cont = ss(A_aug, B_aug, eye(6), zeros(6,2));
sys_aug_disc = c2d(sys_aug_cont, Ts, 'tustin');
Ad = sys_aug_disc.A;
Bd = sys_aug_disc.B;

% Recompute discrete LQR gains
K_aug_d = dlqr(Ad, Bd, Q, R_lq);
Kx_d = K_aug_d(:, 1:4);
Ke_d = K_aug_d(:, 5:6);

disp('Discrete Kx:'); disp(Kx_d)
disp('Discrete Ke:'); disp(Ke_d)

% Save all params
save('LQI_params.mat', 'Kx', 'Ke', 'Kx_d', 'Ke_d', 'x_bar', 'u_bar', 'Ad', 'Bd', 'Ts', 'Q', 'R_lq');