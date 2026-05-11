S = (4.44 / 2)^2 * pi;

v_pump = [];
k1_pump = [];
k2_pump = [];

files = dir('../Data/');
pattern = '^\d+_Pump_BottomTanks_\d+\.mat$';
matches = ~cellfun(@isempty, regexp({files.name}, pattern, 'match'));
pump_files = files(matches);
n_points = numel(pump_files);

for ii = 1:n_points
    load(sprintf('../Data/%s', pump_files(ii).name));
    vv = max(data(2, :));

    time = decimate_avg(data(1, :), 25, 2);
    V1 = decimate_avg(data(2, :), 25, 2);
    V2 = decimate_avg(data(3, :), 25, 2);
    % h1 = decimate_avg(data(4, :), 25, 2);
    h2 = decimate_avg(data(5, :), 25, 2);
    % h3 = decimate_avg(data(6, :), 25, 2);
    h4 = decimate_avg(data(7, :), 25, 2);
    EnFlag = decimate_avg(data(8, :), 25, 2);
    Ts = time(2) - time(1);
    
    idx0 = find(EnFlag > 0, 1, 'first');
    idx0r_2 = idx0 + find(h2(idx0:end) - h2(idx0) > 0.075, 1, 'first');
    idx0r_4 = idx0 + find(h4(idx0:end) - h4(idx0) > 0.075, 1, 'first');
    idx0r_all = [idx0r_2, idx0r_4];
    
    idx0r = max(idx0r_all);
    idx1 = idx0r + find(EnFlag(idx0r:end) == 0, 1, 'first');
    T0 = time(idx0);
    T0r = time(idx0r);
    T1 = time(idx1);

    h2_0 = mean(h2(idx0r-5:idx0r-1));
    h4_0 = mean(h4(idx0r-5:idx0r-1));
    h2_T = mean(h2(idx1:idx1+5));
    h4_T = mean(h4(idx1:idx1+5));
    dV2 = S * (h2_T - h2_0);
    dV4 = S * (h4_T - h4_0);
    F1 = vv * (T1 - T0r);
    F2 = vv * (T1 - T0r);

    theta_2 = dV2 / F1;
    theta_4 = dV4 / F2;

    kappa_2 = theta_4;
    kappa_1 = theta_2;
    
    v_pump = [ v_pump, vv ];
    k1_pump = [ k1_pump, kappa_1 ];
    k2_pump = [ k2_pump, kappa_2 ];
end

[v_pump_sorted, sorting] = sort(v_pump);
[v_pump, ~, groups] = unique(v_pump_sorted);
k1_pump = aggregate_groups(groups, k1_pump(:, sorting));
k2_pump = aggregate_groups(groups, k2_pump(:, sorting));

figure;
plot(v_pump, k1_pump, 'b-o');
hold on;
plot(v_pump, k2_pump, 'r-*');
xlabel('Voltage [V]');
title('Pump gain');
legend('Pump 1', 'Pump 2');
xlim([4, 18]);
ylim([0, 5]);