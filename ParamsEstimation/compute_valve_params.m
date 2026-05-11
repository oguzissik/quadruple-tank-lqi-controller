addpath utils;
idxoff = 1;

S = (4.44 / 2)^2 * pi;
g = 981;
Dsm = 0.31750;
Dmd = 0.47625;
Dlg = 0.55563;
tanks_ids = [];
AA_int = [];
AA_fit = [];
AA_id = [];
DD_norm = [ Dsm, Dlg, Dsm, Dlg ];
AA_nom = (DD_norm ./ 2).^2 * pi;

files = dir('../Data/');
pattern = '^\d+_TankDischarge_\d+\.mat$';
matches = ~cellfun(@isempty, regexp({files.name}, pattern, 'match'));
disch_files = files(matches);
n_points = numel(disch_files);

for ii=1:n_points
    load(sprintf('../Data/%s', disch_files(ii).name));
    time = decimate_avg(data(1, :), 25, 2);
    [ ~, id ] = max(data(2:5, 1));
    
    tanks_ids = [ tanks_ids,  id ];
    hi = decimate_avg(data(id+1, :), 25, 2);
    [xData, yData] = prepareCurveData( time, diff(log([ hi(1), hi])) );

    ft = fittype( 'smoothingspline' );
    opts = fitoptions( 'Method', 'SmoothingSpline' );
    opts.SmoothingParam = 0.99;
    [fitresult, ~] = fit( xData, yData, ft, opts );
    yfilt = fitresult(xData);
    idx0 = idxoff + find(yfilt(idxoff:end) <= -5e-4, 1, 'first');
    idx1 = idx0 + find(yfilt(idx0:end) >= -5e-4, 1, 'first');
    if isempty(idx1)
        idx1 = numel(yfilt);
    end

    T0 = time(idx0);
    T1 = time(idx1);
    Ts = time(idx0+1) - time(idx0);
    
    % Estimate via integration
    hi_avg = (hi(idx0:idx1-1) + hi(idx0+1:idx1)) ./ 2;
    dVi = S * (hi(idx1) - hi(idx0));
    Ai_int = -dVi ./ (sqrt(2*g) * sum(real(sqrt(hi_avg))) * Ts);
    AA_int = [AA_int, Ai_int];

    % Estimate via curve fitting
    hi_ts = max(hi(1, idx0:idx1).', 0*hi(1, idx0:idx1).');
    hi_tilde = sqrt(hi_ts) ;
    t_ts = (0:Ts:Ts*(numel(hi_tilde)-1)).';
    X = sqrt(2*g) / (2*S) .* t_ts;
    Y = hi_tilde(1) - hi_tilde;
    Ai_fit = lsqr(X, Y);
    AA_fit = [ AA_fit, Ai_fit ];

    % Estimate via ID toolbox
%     hi_iddata = iddata(hi_ts, [], Ts, 'Name', 'Level');
%     nlsys = idnlgrey('tank_dynamics', [1, 0, 1], {AA(ii)}, hi_ts(1), 0);
%     setpar(nlsys, 'Fixed', {false});
%     opt = nlgreyestOptions('Display','on');
%     nlsys2 = nlgreyest(hi_iddata, nlsys, opt);
%     Ai_id = nlsys2.Parameters(1).Value;
%     AA_id = [ AA_id, Ai_id ];
end

[ids_sorted, sorting] = sort(tanks_ids);
[tank_ids, ~, groups] = unique(ids_sorted);
AA_int = aggregate_groups(groups, AA_int(:, sorting));
AA_fit = aggregate_groups(groups, AA_fit(:, sorting));

DD_int = 2 .* sqrt(AA_int / pi);
DD_fit = 2 .* sqrt(AA_fit / pi);
% DD_int = 2 .* sqrt(AA_id / pi);

figure;
bar([ DD_norm.', DD_int.', DD_fit.']);
bar([ DD_norm.', DD_int.']);
legend({'Nominal', 'Integration', 'Fitting'}, 'Location','northwest');
xlabel('Tank');
ylabel('cm');
set(gca, 'FontSize', 17);
title('Estimated valve diameter');


figure;
bar([ AA_nom.', AA_int.']);
legend({'Nominal', 'Integration'}, 'Location','northwest');
xlabel('Tank');
ylabel('cm^2');
set(gca, 'FontSize', 17);
title('Estimated cross-section');

diameter_change = 100 * ( DD_fit ./ DD_norm - 1 ) %#ok<NOPTS> 
area_achange = 100 * ( AA_int ./ AA_nom - 1 ) %#ok<NOPTS> 

% Overwrite AA and DD
AA = AA_int;
DD = DD_int;
save('valve_opening', 'AA', 'DD');
