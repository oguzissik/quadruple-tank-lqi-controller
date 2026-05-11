addpath utils;

S = (4.44 / 2)^2 * pi;

tpr_pump = [];
v_pump = [];
k1_pump = [];
k2_pump = [];
gamma1_pump = [];
gamma2_pump = [];

files = dir('../Data/');
pattern = '^\d+_Pump_\d+\.mat$';
matches = ~cellfun(@isempty, regexp({files.name}, pattern, 'match'));
pump_files = files(matches);
n_points = numel(pump_files);

for ii = 1:n_points
    load(sprintf('../Data/%s', pump_files(ii).name));
    vv = max(data(2, :));

    time = decimate_avg(data(1, :), 25, 2);
    V1 = decimate_avg(data(2, :), 25, 2);
    V2 = decimate_avg(data(3, :), 25, 2);
    h1 = decimate_avg(data(4, :), 25, 2);
    h2 = decimate_avg(data(5, :), 25, 2);
    h3 = decimate_avg(data(6, :), 25, 2);
    h4 = decimate_avg(data(7, :), 25, 2);
    EnFlag = decimate_avg(data(8, :), 25, 2);
    Ts = time(2) - time(1);
    
    idx0 = find(EnFlag > 0, 1, 'first');
    idx0r_1 = idx0 + find(h1(idx0:end) - h1(idx0) > 0.075, 1, 'first');
    idx0r_2 = idx0 + find(h2(idx0:end) - h2(idx0) > 0.075, 1, 'first');
    idx0r_3 = idx0 + find(h3(idx0:end) - h3(idx0) > 0.075, 1, 'first');
    idx0r_4 = idx0 + find(h4(idx0:end) - h4(idx0) > 0.075, 1, 'first');
    idx0r_all = [idx0r_1, idx0r_2, idx0r_3, idx0r_4];
    
    idx0r = max(idx0r_all);
    idx1 = idx0r + find(EnFlag(idx0r:end) == 0, 1, 'first');
    T0 = time(idx0);
    T0r = time(idx0r);
    T1 = time(idx1);

    h1_0 = mean(h1(idx0r-5:idx0r-1));
    h2_0 = mean(h2(idx0r-5:idx0r-1));
    h3_0 = mean(h3(idx0r-5:idx0r-1));
    h4_0 = mean(h4(idx0r-5:idx0r-1));
    h1_T = mean(h1(idx1:idx1+5));
    h2_T = mean(h2(idx1:idx1+5));
    h3_T = mean(h3(idx1:idx1+5));
    h4_T = mean(h4(idx1:idx1+5));
    dV1 = S * (h1_T - h1_0);
    dV2 = S * (h2_T - h2_0);
    dV3 = S * (h3_T - h3_0);
    dV4 = S * (h4_T - h4_0);
    F1 = vv * (T1 - T0r);
    F2 = vv * (T1 - T0r);

    theta_1 = dV1 / F2;
    theta_2 = dV2 / F1;
    theta_3 = dV3 / F1;
    theta_4 = dV4 / F2;

    kappa_2 = theta_4 + theta_1;
    gamma_2 = theta_1 / kappa_2;
    kappa_1 = theta_2 + theta_3;
    gamma_1 = theta_3 / kappa_1;
    
    tt = NaN(4, 1);

    for ii=1:4
        idxi = eval(sprintf('idx0r_%d', ii));
        if ~isempty(idxi)
            tt(ii) = time(idxi) - time(idx0);
        end
    end

    v_pump = [ v_pump, vv ];
    k1_pump = [ k1_pump, kappa_1 ];
    k2_pump = [ k2_pump, kappa_2 ];
    gamma1_pump = [ gamma1_pump, gamma_1 ];
    gamma2_pump = [ gamma2_pump, gamma_2 ];
    tpr_pump = [ tpr_pump, tt ];
end

[v_pump_sorted, sorting] = sort(v_pump);
[v_pump, ~, groups] = unique(v_pump_sorted);
k1_pump = aggregate_groups(groups, k1_pump(:, sorting));
k2_pump = aggregate_groups(groups, k2_pump(:, sorting));
gamma1_pump = aggregate_groups(groups, gamma1_pump(:, sorting));
gamma2_pump = aggregate_groups(groups, gamma2_pump(:, sorting));
tpr_pump = aggregate_groups(groups, tpr_pump(:, sorting));
tpr_pump = fill_nans(tpr_pump, 2);

% Fitting via clamped splines
k1_fcn = get_clamped_approximant(v_pump, k1_pump, [4, max(v_pump)], [0, nan], 'spline');
k2_fcn = get_clamped_approximant(v_pump, k2_pump, [4, max(v_pump)], [0, nan], 'spline');
gamma1_fcn = get_clamped_approximant(v_pump, gamma1_pump, [4, max(v_pump)], [0, nan], 'sin');
gamma2_fcn = get_clamped_approximant(v_pump, gamma2_pump, [4, max(v_pump)], [0, nan], 'sin');

% xFit = [ (1-gamma2_pump) .* k2_pump .* v_pump; ...
%          gamma1_pump .* k1_pump .* v_pump; ...
%          (1-gamma1_pump) .* k1_pump .* v_pump; ...
%          gamma2_pump .* k2_pump .* v_pump ];

% A1 = 1 ./ ((1-gamma2_pump) .* k2_pump .* v_pump).';
% A2 = 1 ./ (gamma1_pump .* k1_pump .* v_pump).';
% A3 = 1 ./ ((1-gamma1_pump) .* k1_pump .* v_pump).';
% A4 = 1 ./ (gamma2_pump .* k2_pump .* v_pump).';
% b1 = tpr_pump(1, :).';
% b2 = tpr_pump(2, :).';
% b3 = tpr_pump(3, :).';
% b4 = tpr_pump(4, :).';
% 
% A1 = A1(3:end, :);
% A2 = A2(3:end, :);
% A3 = A3(3:end, :);
% A4 = A4(3:end, :);
% b1 = b1(3:end, :);
% b2 = b2(3:end, :);
% b3 = b3(3:end, :);
% b4 = b4(3:end, :);
% 
% m1 = fitlm(A1, b1);
% m2 = fitlm(A2, b2);
% m3 = fitlm(A3, b3);
% m4 = fitlm(A4, b4);

figure;
v_cont = 0:0.01:18;
plot([0, 18], [0.36, 0.36], 'k--', 'LineWidth', 1.5);
hold on;
plot(v_pump, gamma1_pump, 'b-o', 'LineWidth', 2);
plot(v_pump, gamma2_pump, 'r-*', 'LineWidth', 2);
plot(v_cont, gamma1_fcn(v_cont), 'b:', 'LineWidth', 2);
plot(v_cont, gamma2_fcn(v_cont), 'r:', 'LineWidth', 2);
xlabel('Voltage [V]');
title('Split coefficient');
legend('Nominal', 'Pump 1', 'Pump 2');
set(gca, 'FontSize', 17);
xlim([3, 18]);
ylim([0, 1]);

figure;
plot([0, 18], [3.3, 3.3], 'k--', 'LineWidth', 1.4);
hold on;
plot(v_pump, k1_pump, 'b-o', 'LineWidth', 2);
plot(v_pump, k2_pump, 'r-*', 'LineWidth', 2);
plot(v_cont, k1_fcn(v_cont), 'b:', 'LineWidth', 2);
plot(v_cont, k2_fcn(v_cont), 'r:', 'LineWidth', 2);
xlabel('Voltage [V]');
title('Pump gain');
legend('Nominal', 'Pump 1', 'Pump 2');
set(gca, 'FontSize', 17);
xlim([3, 18]);
ylim([0, 5]);

save('pump_params','v_pump', 'k1_pump', 'k2_pump', 'gamma1_pump', 'gamma2_pump', 'tpr_pump', ...
     'k1_fcn', 'k2_fcn', 'gamma1_fcn', 'gamma2_fcn');