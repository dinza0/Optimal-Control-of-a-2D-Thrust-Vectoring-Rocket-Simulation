%% Part 5.3-5.4: Finite-Horizon LQR for Hover and Weak Descent Linearizations
% Optimal Control of Variable-Mass Thrust-Vectoring Lander
%
% This script:
% 1) Defines hover and weak-descent linearized models
% 2) Solves the finite-horizon Riccati differential equation backward
% 3) Simulates the optimal closed-loop dynamics forward
% 4) Computes state/control/terminal costs
% 5) Checks actuator/state constraints
% 6) Validates PMP stationarity and terminal costate conditions
% 7) Generates state, control, trajectory, Hamiltonian, and comparison plots

clear; clc; close all;

%% =========================
%  1. User Configuration
%  =========================

cfg.g      = 9.81;       % gravity [m/s^2]
cfg.I      = 10;         % moment of inertia [kg*m^2]
cfg.m0     = 20;         % nominal initial mass [kg]
cfg.alpha  = 0.005;      % fuel consumption coefficient [kg/(N*s)]

cfg.Tmin   = 0;          % minimum thrust [N]
cfg.Tmax   = 300;        % maximum thrust [N]
cfg.taumin = -30;        % minimum torque [N*m]
cfg.taumax = 30;         % maximum torque [N*m]
cfg.mdry   = 5;          % dry mass lower bound [kg]

cfg.tf     = 5.0;        % fixed final time [s]
cfg.n      = 7;          % number of states
cfg.r      = 2;          % number of controls

% Nominal hover reference
cfg.hover_px0 = 0;
cfg.hover_pz0 = 50;      % hover altitude [m]
cfg.m_hover   = cfg.m0;

% Nominal weak-descent reference
cfg.weak_px0 = 0;
cfg.weak_pz0 = 50;
cfg.weak_vx  = 10;       % large nominal horizontal velocity [m/s]
cfg.weak_vz  = -3;       % descending nominal vertical velocity [m/s]

% Initial perturbation state:
% x = [dpx; dpz; dvx; dvz; dtheta; domega; dm]
x0 = [ ...
    8;                 % horizontal position error [m]
   -4;                 % vertical position error [m]
    2;                 % horizontal velocity error [m/s]
   -1;                 % vertical velocity error [m/s]
    deg2rad(5);         % pitch angle error [rad]
    deg2rad(2);         % angular velocity error [rad/s]
    0];                % mass perturbation [kg]

% Time grid for plotting and post-processing
cfg.Nplot = 1001;
tgrid = linspace(0, cfg.tf, cfg.Nplot);

%% =========================
%  2. Weight Matrices
%  =========================
% Starting weights from our derivation/discussion.
% You can tune these later for Section 5.5.

Q  = diag([10, 10, 1, 1, 50, 5, 0.1]);
R  = diag([0.05, 2.7]);
Qf = diag([100, 100, 10, 10, 100, 20, 1]);

% Basic checks
assert(all(eig(Q)  >= -1e-12), 'Q must be positive semidefinite.');
assert(all(eig(Qf) >= -1e-12), 'Qf must be positive semidefinite.');
assert(all(eig(R)  >  0),      'R must be positive definite.');

%% =========================
%  3. Define Models
%  =========================

% Hover model: frozen-mass linearization
hover.name = 'Hover';
hover.Afun = @(t) makeA_hover(cfg.g, cfg.m_hover);
hover.Bfun = @(t) makeB_hover(cfg.m_hover, cfg.I, cfg.alpha);
hover.ref  = @(t) reference_hover(t, cfg);

% Weak-descent LTV model
weak.name = 'Weak LTV Descent';
weak.Afun = @(t) makeA_hover(cfg.g, mstar_weak(t, cfg));
weak.Bfun = @(t) makeB_hover(mstar_weak(t, cfg), cfg.I, cfg.alpha);
weak.ref  = @(t) reference_weak(t, cfg);

models = {hover, weak};

%% =========================
%  4. Solve LQR for Both Models
%  =========================

results = cell(numel(models), 1);

for i = 1:numel(models)
    model = models{i};

    fprintf('\n======================================================\n');
    fprintf('Solving finite-horizon LQR for: %s\n', model.name);
    fprintf('======================================================\n');

    results{i} = solve_lqr_model(model, cfg, Q, R, Qf, x0, tgrid);

    print_summary(results{i}, cfg);
end

%% =========================
%  5.5 Weight Tuning Sweep
%  =========================

weightCases = {};

weightCases{end+1}.name = 'Aggressive infeasible';
weightCases{end}.Q  = diag([10, 10, 1, 1, 100, 10, 0.1]);
weightCases{end}.R  = diag([0.01, 0.01]);
weightCases{end}.Qf = diag([100, 100, 10, 10, 500, 50, 1]);

weightCases{end+1}.name = 'Moderate torque penalty';
weightCases{end}.Q  = diag([10, 10, 1, 1, 50, 5, 0.1]);
weightCases{end}.R  = diag([0.05, 1.0]);
weightCases{end}.Qf = diag([100, 100, 10, 10, 100, 20, 1]);

weightCases{end+1}.name = 'Tuned feasible';
weightCases{end}.Q  = diag([10, 10, 1, 1, 50, 5, 0.1]);
weightCases{end}.R  = diag([0.05, 2.7]);
weightCases{end}.Qf = diag([100, 100, 10, 10, 100, 20, 1]);

weightCases{end+1}.name = 'Conservative';
weightCases{end}.Q  = diag([10, 10, 1, 1, 50, 5, 0.1]);
weightCases{end}.R  = diag([0.10, 5.0]);
weightCases{end}.Qf = diag([100, 100, 10, 10, 100, 20, 1]);

x0_base = [8; -4; 2; -1; deg2rad(5); deg2rad(2); 0];

summaryRows = [];

for c = 1:numel(weightCases)
    Q  = weightCases{c}.Q;
    R  = weightCases{c}.R;
    Qf = weightCases{c}.Qf;

    for i = 1:numel(models)
        res = solve_lqr_model(models{i}, cfg, Q, R, Qf, x0_base, tgrid);

        feasible = ~(res.thrust_violation || res.torque_violation || ...
            res.altitude_violation || res.mass_violation);

        summaryRows = [summaryRows;
            {weightCases{c}.name, res.name, res.Jtotal, res.Jcontrol, ...
            norm(res.x(:,end)), res.max_Tphys, max(abs(res.tauphys)), feasible}];
    end
end

weightTable = cell2table(summaryRows, ...
    'VariableNames', {'Case','Model','J_total','J_control','TerminalNorm', ...
    'MaxThrust','MaxAbsTorque','Feasible'});

disp(weightTable);
writetable(weightTable, 'part55_weight_tuning_table.csv');

%% =========================
%  5.5 Initial Condition Sensitivity Sweep
%  =========================

Q_tuned  = diag([10, 10, 1, 1, 50, 5, 0.1]);
R_tuned  = diag([0.05, 2.7]);
Qf_tuned = diag([100, 100, 10, 10, 100, 20, 1]);

icCases = {};

icCases{end+1}.name = 'Nominal IC';
icCases{end}.x0 = x0_base;

icCases{end+1}.name = '25 percent larger IC';
icCases{end}.x0 = 1.25*x0_base;

icCases{end+1}.name = '50 percent larger IC';
icCases{end}.x0 = 1.50*x0_base;

icCases{end+1}.name = 'High vx error';
icCases{end}.x0 = [8; -4; 5; -1; deg2rad(5); deg2rad(2); 0];

icCases{end+1}.name = 'High theta error';
icCases{end}.x0 = [8; -4; 2; -1; deg2rad(10); deg2rad(2); 0];

summaryRows = [];

for c = 1:numel(icCases)
    x0_test = icCases{c}.x0;

    for i = 1:numel(models)
        res = solve_lqr_model(models{i}, cfg, Q_tuned, R_tuned, Qf_tuned, x0_test, tgrid);

        feasible = ~(res.thrust_violation || res.torque_violation || ...
            res.altitude_violation || res.mass_violation);

        peakStateNorm = max(vecnorm(res.x,2,1));

        summaryRows = [summaryRows;
            {icCases{c}.name, res.name, norm(x0_test), norm(res.x(:,end)), ...
            peakStateNorm, res.Jtotal, res.Jcontrol, res.max_Tphys, ...
            max(abs(res.tauphys)), feasible}];
    end
end

icTable = cell2table(summaryRows, ...
    'VariableNames', {'IC_Case','Model','InitialNorm','TerminalNorm', ...
    'PeakStateNorm','J_total','J_control', ...
    'MaxThrust','MaxAbsTorque','Feasible'});

disp(icTable);
writetable(icTable, 'part55_ic_sensitivity_table.csv');

figure;
hold on; grid on;

for c = 1:numel(icCases)
    x0_test = icCases{c}.x0;
    res = solve_lqr_model(models{1}, cfg, Q_tuned, R_tuned, Qf_tuned, x0_test, tgrid);
    plot(res.t, vecnorm(res.x,2,1), 'LineWidth', 1.5);
end

xlabel('Time [s]');
ylabel('||x(t)||_2');
title('Hover LQR Initial Condition Sensitivity');
legend({icCases{1}.name, icCases{2}.name, icCases{3}.name, ...
    icCases{4}.name, icCases{5}.name}, 'Location', 'best');

figure;
hold on; grid on;

for c = 1:numel(icCases)
    x0_test = icCases{c}.x0;
    res = solve_lqr_model(models{2}, cfg, Q_tuned, R_tuned, Qf_tuned, x0_test, tgrid);
    plot(res.t, vecnorm(res.x,2,1), 'LineWidth', 1.5);
end

xlabel('Time [s]');
ylabel('||x(t)||_2');
title('Weakly Controlled LQR Initial Condition Sensitivity');
legend({icCases{1}.name, icCases{2}.name, icCases{3}.name, ...
    icCases{4}.name, icCases{5}.name}, 'Location', 'best');

%% =========================
%  5. Comparison Table
%  =========================

fprintf('\n\n================ MODEL COMPARISON ================\n');
fprintf('%-20s  %-12s  %-12s  %-12s  %-12s  %-12s\n', ...
    'Model', 'J_total', 'J_state', 'J_control', 'J_final', '||x(tf)||');

for i = 1:numel(results)
    res = results{i};
    fprintf('%-20s  %-12.4e  %-12.4e  %-12.4e  %-12.4e  %-12.4e\n', ...
        res.name, res.Jtotal, res.Jstate, res.Jcontrol, res.Jfinal, norm(res.x(:,end)));
end

%% =========================
%  6. Plot Results
%  =========================

plot_results(results, cfg);

%% ============================================================
%  Local Functions
%  ============================================================

function res = solve_lqr_model(model, cfg, Q, R, Qf, x0, tgrid)

    n = cfg.n;

    % -------------------------
    % Solve Riccati equation backward
    % -------------------------
    Ptf = Qf;
    p0_vec = Ptf(:);

    riccati_rhs = @(t, pvec) riccati_ode(t, pvec, model.Afun, model.Bfun, Q, R, n);

    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
    [tP_raw, Pvec_raw] = ode45(riccati_rhs, [cfg.tf 0], p0_vec, opts);

    % Sort Riccati solution in ascending time for interpolation
    [tP, sortIdx] = sort(tP_raw);
    Pvec = Pvec_raw(sortIdx, :);

    % Evaluate P(t) on plotting grid
    Pgrid = zeros(n, n, numel(tgrid));
    for k = 1:numel(tgrid)
        Pgrid(:,:,k) = getP(tgrid(k), tP, Pvec, n);
    end

    % -------------------------
    % Simulate closed-loop system forward
    % -------------------------
    closed_loop_rhs = @(t, x) closed_loop_ode(t, x, model.Afun, model.Bfun, R, tP, Pvec, n);

    [tx_raw, x_raw] = ode45(closed_loop_rhs, [0 cfg.tf], x0, opts);

    % Interpolate state to uniform grid
    xgrid = interp1(tx_raw, x_raw, tgrid, 'pchip').';
    
    % -------------------------
    % Compute controls, costates, Hamiltonian, references
    % -------------------------
    ugrid      = zeros(cfg.r, numel(tgrid));
    lambdagrid = zeros(n, numel(tgrid));
    Hgrid      = zeros(1, numel(tgrid));
    stationarity_res = zeros(cfg.r, numel(tgrid));

    xabs_grid  = zeros(n, numel(tgrid));
    Tphys      = zeros(1, numel(tgrid));
    tauphys    = zeros(1, numel(tgrid));

    for k = 1:numel(tgrid)
        t = tgrid(k);
        x = xgrid(:,k);

        A = model.Afun(t);
        B = model.Bfun(t);
        P = Pgrid(:,:,k);

        K = R \ (B.' * P);
        u = -K * x;

        lambda = 2 * P * x;

        xdot = A*x + B*u;

        H = x.'*Q*x + u.'*R*u + lambda.'*xdot;

        ru = 2*R*u + B.'*lambda;

        ref = model.ref(t);
        xabs = ref.x + x;

        if strcmp(model.name, 'Hover')
            Tnom = cfg.m_hover * cfg.g;
            taunom = 0;
        else
            Tnom = ref.Tnom;
            taunom = 0;
        end

        Tphys(k)   = Tnom + u(1);
        tauphys(k) = taunom + u(2);

        ugrid(:,k) = u;
        lambdagrid(:,k) = lambda;
        Hgrid(k) = H;
        stationarity_res(:,k) = ru;
        xabs_grid(:,k) = xabs;
    end

    % -------------------------
    % Cost calculation
    % -------------------------
    state_integrand   = zeros(1, numel(tgrid));
    control_integrand = zeros(1, numel(tgrid));

    for k = 1:numel(tgrid)
        x = xgrid(:,k);
        u = ugrid(:,k);

        state_integrand(k)   = x.'*Q*x;
        control_integrand(k) = u.'*R*u;
    end

    Jstate   = trapz(tgrid, state_integrand);
    Jcontrol = trapz(tgrid, control_integrand);
    Jfinal   = xgrid(:,end).'*Qf*xgrid(:,end);
    Jtotal   = Jstate + Jcontrol + Jfinal;

    % -------------------------
    % Validation checks
    % -------------------------
    lambda_tf_expected = 2*Qf*xgrid(:,end);
    lambda_tf_actual   = lambdagrid(:,end);
    terminal_costate_residual = lambda_tf_actual - lambda_tf_expected;

    symmetry_error = zeros(1, numel(tgrid));
    for k = 1:numel(tgrid)
        P = Pgrid(:,:,k);
        symmetry_error(k) = norm(P - P.', 'fro');
    end

    % Actuator and state constraints
    thrust_violation_low  = min(Tphys) < cfg.Tmin;
    thrust_violation_high = max(Tphys) > cfg.Tmax;
    torque_violation_low  = min(tauphys) < cfg.taumin;
    torque_violation_high = max(tauphys) > cfg.taumax;

    pz_min = min(xabs_grid(2,:));
    m_min  = min(xabs_grid(7,:));

    altitude_violation = pz_min < 0;
    mass_violation     = m_min < cfg.mdry;

    % Pack results
    res.name = model.name;
    res.t = tgrid;
    res.x = xgrid;
    res.u = ugrid;
    res.lambda = lambdagrid;
    res.P = Pgrid;
    res.H = Hgrid;
    res.xabs = xabs_grid;
    res.Tphys = Tphys;
    res.tauphys = tauphys;

    res.Jstate = Jstate;
    res.Jcontrol = Jcontrol;
    res.Jfinal = Jfinal;
    res.Jtotal = Jtotal;

    res.state_integrand = state_integrand;
    res.control_integrand = control_integrand;

    res.stationarity_res = stationarity_res;
    res.terminal_costate_residual = terminal_costate_residual;
    res.symmetry_error = symmetry_error;

    res.max_stationarity_res_norm = max(vecnorm(stationarity_res, 2, 1));
    res.terminal_costate_res_norm = norm(terminal_costate_residual);
    res.max_symmetry_error = max(symmetry_error);

    res.max_deltaT = max(abs(ugrid(1,:)));
    res.max_tau    = max(abs(ugrid(2,:)));
    res.energy_uTu = trapz(tgrid, sum(ugrid.^2, 1));

    res.min_Tphys = min(Tphys);
    res.max_Tphys = max(Tphys);
    res.min_tauphys = min(tauphys);
    res.max_tauphys = max(tauphys);
    res.min_pz = pz_min;
    res.min_m  = m_min;

    res.thrust_violation = thrust_violation_low || thrust_violation_high;
    res.torque_violation = torque_violation_low || torque_violation_high;
    res.altitude_violation = altitude_violation;
    res.mass_violation = mass_violation;
end

function dpvec = riccati_ode(t, pvec, Afun, Bfun, Q, R, n)
    P = reshape(pvec, n, n);
    P = 0.5*(P + P.');  % numerical symmetry cleanup

    A = Afun(t);
    B = Bfun(t);

    % Riccati equation:
    % -Pdot = A'P + PA - PBR^-1B'P + Q
    Pdot = -(A.'*P + P*A - P*B*(R \ B.')*P + Q);

    dpvec = Pdot(:);
end

function dx = closed_loop_ode(t, x, Afun, Bfun, R, tP, Pvec, n)
    A = Afun(t);
    B = Bfun(t);
    P = getP(t, tP, Pvec, n);

    K = R \ (B.' * P);
    u = -K*x;

    dx = A*x + B*u;
end

function P = getP(t, tP, Pvec, n)
    p_interp = interp1(tP, Pvec, t, 'pchip');
    P = reshape(p_interp.', n, n);
    P = 0.5*(P + P.');
end

function A = makeA_hover(g, m)
    A = [0 0 1 0 0 0 0;
         0 0 0 1 0 0 0;
         0 0 0 0 g 0 0;
         0 0 0 0 0 0 -g/m;
         0 0 0 0 0 1 0;
         0 0 0 0 0 0 0;
         0 0 0 0 0 0 0];
end

function B = makeB_hover(m, I, alpha)
    B = [0 0;
         0 0;
         0 0;
         1/m 0;
         0 0;
         0 1/I;
         -alpha 0];
end

function m = mstar_weak(t, cfg)
    m = cfg.m0 * exp(-cfg.alpha * cfg.g * t);
end

function ref = reference_hover(t, cfg)
    ref.x = [cfg.hover_px0;
             cfg.hover_pz0;
             0;
             0;
             0;
             0;
             cfg.m_hover];

    ref.Tnom = cfg.m_hover * cfg.g;
end

function ref = reference_weak(t, cfg)
    m = mstar_weak(t, cfg);

    ref.x = [cfg.weak_px0 + cfg.weak_vx*t;
             cfg.weak_pz0 + cfg.weak_vz*t;
             cfg.weak_vx;
             cfg.weak_vz;
             0;
             0;
             m];

    ref.Tnom = m * cfg.g;
end

function print_summary(res, cfg)

    fprintf('\n--- Cost Results: %s ---\n', res.name);
    fprintf('J_state                  = %.6e\n', res.Jstate);
    fprintf('J_control                = %.6e\n', res.Jcontrol);
    fprintf('J_final                  = %.6e\n', res.Jfinal);
    fprintf('J_total                  = %.6e\n', res.Jtotal);
    fprintf('||x(tf)||                = %.6e\n', norm(res.x(:,end)));
    fprintf('Integral u^T u dt         = %.6e\n', res.energy_uTu);

    fprintf('\n--- Control Usage ---\n');
    fprintf('max |delta T|             = %.6e N\n', res.max_deltaT);
    fprintf('max |tau|                 = %.6e N*m\n', res.max_tau);
    fprintf('min physical T            = %.6e N\n', res.min_Tphys);
    fprintf('max physical T            = %.6e N\n', res.max_Tphys);
    fprintf('min physical tau          = %.6e N*m\n', res.min_tauphys);
    fprintf('max physical tau          = %.6e N*m\n', res.max_tauphys);

    fprintf('\n--- Physical Constraints ---\n');
    fprintf('min altitude pz           = %.6e m\n', res.min_pz);
    fprintf('min mass m                = %.6e kg\n', res.min_m);

    if res.thrust_violation
        fprintf('WARNING: thrust bound violated. Bounds: [%.2f, %.2f] N\n', cfg.Tmin, cfg.Tmax);
    else
        fprintf('Thrust bounds satisfied.\n');
    end

    if res.torque_violation
        fprintf('WARNING: torque bound violated. Bounds: [%.2f, %.2f] N*m\n', cfg.taumin, cfg.taumax);
    else
        fprintf('Torque bounds satisfied.\n');
    end

    if res.altitude_violation
        fprintf('WARNING: altitude constraint pz >= 0 violated.\n');
    else
        fprintf('Altitude constraint satisfied.\n');
    end

    if res.mass_violation
        fprintf('WARNING: dry mass constraint m >= mdry violated. mdry = %.2f kg\n', cfg.mdry);
    else
        fprintf('Dry mass constraint satisfied.\n');
    end

    fprintf('\n--- PMP / Riccati Validation ---\n');
    fprintf('max ||2Ru + B^T lambda||  = %.6e\n', res.max_stationarity_res_norm);
    fprintf('||lambda(tf)-2Qf x(tf)||  = %.6e\n', res.terminal_costate_res_norm);
    fprintf('max ||P-P^T||_F           = %.6e\n', res.max_symmetry_error);
    fprintf('Hamiltonian initial       = %.6e\n', res.H(1));
    fprintf('Hamiltonian final         = %.6e\n', res.H(end));
    fprintf('Hamiltonian variation     = %.6e\n', max(res.H)-min(res.H));
end

function plot_results(results, cfg)

    labels = {'\delta p_x', '\delta p_z', '\delta v_x', '\delta v_z', '\delta\theta', '\delta\omega', '\delta m'};
    units  = {'m', 'm', 'm/s', 'm/s', 'rad', 'rad/s', 'kg'};

    for i = 1:numel(results)
        res = results{i};
        t = res.t;

        %% Perturbation states
        figure('Name', [res.name ' Perturbation States']);
        tiledlayout(4,2);

        for k = 1:7
            nexttile;
            plot(t, res.x(k,:), 'LineWidth', 1.5);
            grid on;
            xlabel('Time [s]');
            ylabel([labels{k} ' [' units{k} ']']);
            title(labels{k});
        end

        %% Physical states
        figure('Name', [res.name ' Physical States']);
        tiledlayout(4,2);

        physical_labels = {'p_x', 'p_z', 'v_x', 'v_z', '\theta', '\omega', 'm'};

        for k = 1:7
            nexttile;
            plot(t, res.xabs(k,:), 'LineWidth', 1.5);
            grid on;
            xlabel('Time [s]');
            ylabel(physical_labels{k});
            title(['Physical ' physical_labels{k}]);
        end

        %% Controls
        figure('Name', [res.name ' Controls']);
        tiledlayout(3,1);

        nexttile;
        plot(t, res.u(1,:), 'LineWidth', 1.5);
        grid on;
        xlabel('Time [s]');
        ylabel('\delta T or \tilde T [N]');
        title('Optimal Thrust Perturbation');

        nexttile;
        plot(t, res.u(2,:), 'LineWidth', 1.5);
        grid on;
        xlabel('Time [s]');
        ylabel('\tau [N*m]');
        title('Optimal Torque');

        nexttile;
        plot(t, res.Tphys, 'LineWidth', 1.5); hold on;
        yline(cfg.Tmin, '--');
        yline(cfg.Tmax, '--');
        grid on;
        xlabel('Time [s]');
        ylabel('T_{phys} [N]');
        title('Physical Thrust with Bounds');

        %% Hamiltonian
        figure('Name', [res.name ' Hamiltonian']);
        plot(t, res.H, 'LineWidth', 1.5);
        grid on;
        xlabel('Time [s]');
        ylabel('Hamiltonian');
        title([res.name ' Hamiltonian Evolution']);

        %% Mission trajectory px vs pz with attitude vectors
        figure('Name', [res.name ' Mission Trajectory']);
        plot(res.xabs(1,:), res.xabs(2,:), 'LineWidth', 1.8); hold on;
        grid on;
        xlabel('p_x [m]');
        ylabel('p_z [m]');
        title([res.name ' Mission Trajectory with Attitude Vectors']);
        axis equal;

        stride = max(1, floor(numel(t)/20));
        scale = 2.0;

        for k = 1:stride:numel(t)
            px = res.xabs(1,k);
            pz = res.xabs(2,k);
            th = res.xabs(5,k);

            dx = scale*sin(th);
            dz = scale*cos(th);

            quiver(px, pz, dx, dz, 0, 'LineWidth', 1.0, 'MaxHeadSize', 1.5);
        end
    end

    %% Cross-model comparison plots
    if numel(results) >= 2
        figure('Name', 'Comparison: Terminal Error and Controls');
        tiledlayout(3,1);

        nexttile;
        hold on; grid on;
        for i = 1:numel(results)
            plot(results{i}.t, vecnorm(results{i}.x,2,1), 'LineWidth', 1.5);
        end
        xlabel('Time [s]');
        ylabel('||x(t)||_2');
        legend(results{1}.name, results{2}.name);
        title('Perturbation State Norm');

        nexttile;
        hold on; grid on;
        for i = 1:numel(results)
            plot(results{i}.t, results{i}.u(1,:), 'LineWidth', 1.5);
        end
        xlabel('Time [s]');
        ylabel('Thrust perturbation [N]');
        legend(results{1}.name, results{2}.name);
        title('Optimal Thrust Perturbation Comparison');

        nexttile;
        hold on; grid on;
        for i = 1:numel(results)
            plot(results{i}.t, results{i}.u(2,:), 'LineWidth', 1.5);
        end
        xlabel('Time [s]');
        ylabel('Torque [N*m]');
        legend(results{1}.name, results{2}.name);
        title('Optimal Torque Comparison');

        figure('Name', 'Comparison: Mission Trajectories');
        hold on; grid on; axis equal;
        for i = 1:numel(results)
            plot(results{i}.xabs(1,:), results{i}.xabs(2,:), 'LineWidth', 1.8);
        end
        xlabel('p_x [m]');
        ylabel('p_z [m]');
        legend(results{1}.name, results{2}.name);
        title('Mission Trajectory Comparison');
    end
end